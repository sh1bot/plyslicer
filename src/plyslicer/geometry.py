from dataclasses import dataclass

import numpy as np
import shapely
import trimesh
from scipy.ndimage import distance_transform_edt


class GeometryError(ValueError):
    pass


class DisconnectedLayerError(GeometryError):
    """A valid solid cannot be laminated in the current direction."""


@dataclass
class Frame:
    origin: np.ndarray
    basis: np.ndarray  # columns: local x, local y, insertion direction

    def local(self, mesh):
        result = mesh.copy()
        result.vertices = (mesh.vertices - self.origin) @ self.basis
        return result

    def world(self, mesh):
        result = mesh.copy()
        result.vertices = mesh.vertices @ self.basis.T + self.origin
        return result


def load_model(path):
    scene = trimesh.load_scene(path, process=True)
    if not scene.geometry:
        raise GeometryError("Input has no mesh geometry.")
    # Apply 3MF scene/build transforms and convert declared units to millimeters.
    if str(path).lower().endswith('.3mf'):
        scene = scene.convert_units('mm')
    mesh = scene.to_mesh()
    if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume <= 0:
        raise GeometryError("Input must be a closed, consistently oriented positive-volume mesh; repair it first.")
    if len(mesh.split()) != 1:
        raise GeometryError("Select one connected solid; multiple solids are not supported.")
    return mesh


def make_frame(mesh, direction):
    direction = np.asarray(direction, dtype=float)
    norm = np.linalg.norm(direction)
    if not np.isfinite(norm) or norm < 1e-12:
        raise GeometryError("Direction must be a finite nonzero vector.")
    z = direction / norm
    seed = np.eye(3)[np.argmin(np.abs(z))]
    x = np.cross(seed, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return Frame(mesh.bounds.mean(axis=0), np.column_stack((x, y, z)))


def projection_check(mesh):
    """Detect repeated entries/exits over positive projected area.

    For a valid solid, more than one exit at a projected point implies an
    undercut. Compare summed face projection area against its planar union.
    This examines triangles, rather than only sampled rays. Floating-point
    tolerance is explicit; this is not an exact-arithmetic mesh validator.
    """
    tri = mesh.triangles[:, :, :2]
    edge1, edge2 = tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]
    cross = edge1[:, 0] * edge2[:, 1] - edge1[:, 1] * edge2[:, 0]
    areas = []
    excess = []
    for sign in (1, -1):
        selected = cross * sign > 0
        polygons = shapely.polygons(tri[selected])
        area = float(shapely.area(shapely.union_all(polygons)))
        total = float(np.abs(cross[selected]).sum() / 2)
        areas.append(area)
        excess.append(max(0., total - area))
    area = min(areas)
    tolerance = max(1e-9, area * 1e-8)
    return {'projected_area_mm2': area, 'overlap_area_mm2': max(excess),
            'area_tolerance_mm2': tolerance,
            'passes': bool(area > 0 and max(excess) <= tolerance)}


def choose_frame(mesh, direction=None):
    if direction is not None:
        candidates = [make_frame(mesh, direction)]
    else:
        # Broad thin shapes have large projected area / mean ray thickness.
        _, vectors = np.linalg.eigh(np.cov(mesh.vertices.T))
        candidates = []
        for d in [*vectors.T, *np.eye(3)]:
            f = make_frame(mesh, d)
            if not any(abs(f.basis[:, 2] @ c.basis[:, 2]) > .99999 for c in candidates):
                candidates.append(f)
        candidates.sort(key=lambda f: float(np.abs(mesh.face_normals @ f.basis[:, 2]) @ mesh.area_faces), reverse=True)
    frame = candidates[0]
    check = projection_check(frame.local(mesh))
    diagnostics = [{'direction': frame.basis[:, 2].tolist(), **check}]
    # Exterior undercuts are allowed: the INTERFACES must be height fields.
    # Ray envelopes can span cavities; downstream connectedness checks reject
    # splits that isolate islands rather than silently exporting fragments.
    return candidates[0], diagnostics


def sample_shell(mesh, pitch, max_points=200000):
    lo, hi = mesh.bounds
    # A surrounding rectangular grid allows closed, nonintersecting cutters.
    nx, ny = np.ceil((hi[:2] - lo[:2]) / pitch).astype(int) + 3
    if nx * ny > max_points:
        raise GeometryError(f"Grid needs {nx * ny} samples; increase --pitch (limit {max_points}).")
    xs = np.linspace(lo[0] - pitch, hi[0] + pitch, nx)
    ys = np.linspace(lo[1] - pitch, hi[1] + pitch, ny)
    xx, yy = np.meshgrid(xs, ys)
    xy = np.column_stack((xx.ravel(), yy.ravel()))
    lower = np.full(len(xy), np.inf)
    upper = np.full(len(xy), -np.inf)
    # Bound ray query memory independently of input and grid size.
    for start in range(0, len(xy), 512):
        points = xy[start:start + 512]
        origins = np.column_stack((points, np.full(len(points), lo[2] - pitch)))
        loc, rays, _ = mesh.ray.intersects_location(origins, np.tile([0., 0., 1.], (len(points), 1)))
        if len(loc) == 0:
            continue
        np.minimum.at(lower, start + rays, loc[:, 2])
        np.maximum.at(upper, start + rays, loc[:, 2])
    valid = np.isfinite(lower) & np.isfinite(upper) & (upper - lower > 1e-8)
    if valid.sum() < 4:
        raise GeometryError("Too few shell samples; decrease --pitch.")
    valid = valid.reshape(ny, nx)
    indices = distance_transform_edt(~valid, return_distances=False, return_indices=True)
    lower = lower.reshape(ny, nx)[tuple(indices)]
    upper = upper.reshape(ny, nx)[tuple(indices)]
    return xy, lower, upper, valid


def slab(xy, bottom, top):
    ny, nx = bottom.shape
    n = nx * ny
    vertices = np.vstack((np.column_stack((xy, bottom.ravel())), np.column_stack((xy, top.ravel()))))
    j, i = np.meshgrid(np.arange(ny - 1), np.arange(nx - 1), indexing='ij')
    a = (j * nx + i).ravel()
    faces = np.vstack((np.column_stack((a, a + 1, a + nx + 1)), np.column_stack((a, a + nx + 1, a + nx))))
    result = [faces[:, ::-1], faces + n]
    edge = np.concatenate((np.arange(nx), np.arange(2 * nx - 1, n, nx),
                           np.arange(n - 2, n - nx - 1, -1), np.arange(n - 2 * nx, 0, -nx)))
    nxt = np.roll(edge, -1)
    result.extend((np.column_stack((edge, nxt, nxt + n)), np.column_stack((edge, nxt + n, edge + n))))
    return trimesh.Trimesh(vertices=vertices, faces=np.vstack(result), process=True)


def split_model(mesh, layers=2, pitch=.8, direction=None, min_thickness=.4,
                pins=0, pin_diameter=2., clearance=.15):
    options = dict(layers=layers, pitch=pitch, min_thickness=min_thickness,
                   pins=pins, pin_diameter=pin_diameter, clearance=clearance)
    try:
        return _split_model(mesh, direction=direction, **options)
    except DisconnectedLayerError as first_error:
        if direction is not None:
            raise
        first_message = str(first_error)
    # A broad projection alone does not imply a connected lamination. Explore
    # tilted insertion directions, keeping every component and all validation.
    # World diagonals help models imported in a different print orientation;
    # PCA diagonals make the search useful independent of that orientation.
    primary, _ = choose_frame(mesh)
    _, axes = np.linalg.eigh(np.cov(mesh.vertices.T))
    directions = [*axes.T, *np.eye(3)]
    for basis in (axes, np.eye(3)):
        for i in range(3):
            for j in range(i + 1, 3):
                for sign in (-1, 1):
                    directions.append(basis[:, i] + sign * basis[:, j])
    candidates = []
    for d in directions:
        frame = make_frame(mesh, d)
        if abs(frame.basis[:, 2] @ primary.basis[:, 2]) > .99999:
            continue
        if not any(abs(frame.basis[:, 2] @ f.basis[:, 2]) > .99999 for f in candidates):
            candidates.append(frame)
    candidates.sort(key=lambda f: float(np.abs(mesh.face_normals @ f.basis[:, 2]) @ mesh.area_faces), reverse=True)
    attempts = [{'direction': primary.basis[:, 2].tolist(), 'error': first_message}]
    for frame in candidates:
        d = frame.basis[:, 2]
        try:
            parts, report = _split_model(mesh, direction=d, **options)
        except GeometryError as exc:
            attempts.append({'direction': d.tolist(), 'error': str(exc)})
            continue
        report['automatic_direction_failures'] = attempts
        report['warnings'].append(f'Used a fallback insertion direction after {len(attempts)} rejected splits.')
        return parts, report
    raise GeometryError('Automatic direction search found no valid connected split. '
                        'Try --direction X,Y,Z or a finer --pitch. '
                        f'First failure: {first_message}')


def _split_model(mesh, layers=2, pitch=.8, direction=None, min_thickness=.4,
                 pins=0, pin_diameter=2., clearance=.15):
    if layers < 2 or layers > 32:
        raise GeometryError("--layers must be between 2 and 32.")
    if not all(np.isfinite(v) and v > 0 for v in [pitch, min_thickness, pin_diameter]):
        raise GeometryError("Pitch, minimum thickness, and pin diameter must be positive finite numbers.")
    if not np.isfinite(clearance) or clearance < 0 or pins < 0 or pins > 8:
        raise GeometryError("Clearance must be nonnegative and --pins between 0 and 8.")
    frame, checks = choose_frame(mesh, direction)
    local = frame.local(mesh)
    xy, lower, upper, valid = sample_shell(local, pitch)
    thickness = (upper - lower)[valid]
    # At rounded perimeter edges thickness tends to zero. Report this honestly;
    # the minimum setting is a bulk (5th-percentile) guard, not an edge guarantee.
    bulk = float(np.percentile(thickness, 5) / layers)
    if bulk < min_thickness:
        raise GeometryError(f"5th-percentile layer thickness {bulk:.3f} mm is below --min-thickness {min_thickness}.")
    lo, hi = local.bounds[:, 2]
    margin = max(hi - lo, pitch) + 1
    surfaces = [np.full(lower.shape, lo - margin)]
    surfaces += [lower + (upper - lower) * k / layers for k in range(1, layers)]
    surfaces += [np.full(lower.shape, hi + margin)]
    parts = []
    for k in range(layers):
        cutter = slab(xy, surfaces[k], surfaces[k + 1])
        part = trimesh.boolean.intersection([local, cutter], engine='manifold')
        if part.is_empty or not part.is_watertight or part.volume <= 0:
            raise GeometryError(f"Layer {k + 1} did not produce a closed solid; try a smaller --pitch.")
        if len(part.split()) != 1:
            raise DisconnectedLayerError(f"Layer {k + 1} is disconnected; try another --direction, fewer layers, or a smaller --pitch.")
        parts.append(part)
    uncut_volume = sum(p.volume for p in parts)
    if abs(uncut_volume - local.volume) > local.volume * 1e-4:
        raise GeometryError("Boolean volume conservation failed; no parts exported.")
    pin_locations = []
    if pins:
        # Select well-separated bore centers with conservative footprint margins.
        ny, nx = valid.shape
        dx = xy[1, 0] - xy[0, 0]
        dy = xy[nx, 1] - xy[0, 1]
        distance = distance_transform_edt(valid, sampling=(dy, dx))
        radius = pin_diameter / 2 + clearance
        eligible = distance > radius + max(2 * pitch, min_thickness)
        separation = np.full(valid.shape, np.inf)
        cylinders = []
        for _ in range(pins):
            score = np.where(eligible, np.minimum(distance, separation), -1.)
            index = int(np.argmax(score))
            if score.ravel()[index] <= radius * 3:
                raise GeometryError("Not enough interior space for the requested alignment bores.")
            point = xy[index]
            pin_locations.append([float(point[0]), float(point[1])])
            cylinder = trimesh.creation.cylinder(radius=radius, height=hi-lo+2*margin, sections=48)
            cylinder.apply_translation([*point, (lo + hi) / 2])
            cylinders.append(cylinder)
            separation = np.minimum(separation, np.linalg.norm(xy-point, axis=1).reshape(valid.shape))
        for k, part in enumerate(parts):
            parts[k] = trimesh.boolean.difference([part, *cylinders], engine='manifold')
            if not parts[k].is_watertight or parts[k].volume <= 0 or len(parts[k].split()) != 1:
                raise GeometryError("Alignment bores broke a layer; reduce their size or count.")
    report = {
        'units': 'mm', 'layers': layers, 'pitch_mm': pitch,
        'frame_origin_world': frame.origin.tolist(), 'frame_basis_columns': frame.basis.tolist(),
        'insertion_direction_world': frame.basis[:, 2].tolist(),
        'direction_checks': checks,
        'assembly': 'Each interface is a single-valued piecewise-linear height field with no undercuts. '
                    'Assemble in increasing layer order, lowering each next layer along the negative insertion direction. '
                    'The original exterior may have overhangs; projection checks describe the input, not the interfaces. '
                    'Numerical booleans and manufacturing tolerances are not a mechanical fit certification.',
        'sampled_axial_thickness_mm': {'min': float(thickness.min()), 'p05': float(np.percentile(thickness, 5)),
                                       'median': float(np.median(thickness)), 'max': float(thickness.max())},
        'bulk_layer_thickness_p05_mm': bulk,
        'warnings': ['Thickness is measured along insertion rays, not surface normals. '
                     'Rounded edges can be thinner than the minimum setting.',
                     'Interfaces approximate sampled shell fractions; reduce pitch to check convergence.',
                     'Curved laminations may need supports. Choose print orientations and inspect slicer previews.'],
        'input_volume_mm3': float(local.volume), 'volume_before_bores_mm3': float(uncut_volume),
        'output_volume_mm3': float(sum(p.volume for p in parts)),
        'alignment': {'bore_centers_local_xy_mm': pin_locations, 'rod_diameter_mm': pin_diameter,
                      'radial_clearance_mm': clearance,
                      'instructions': 'Use removable straight rods through the stack; remove after glue sets. '
                                      'No captive keys or permanent protruding pins are generated.'},
        'parts': [{'layer': i+1, 'volume_mm3': float(p.volume), 'watertight': bool(p.is_watertight)}
                  for i, p in enumerate(parts)]}
    if not checks[0]['passes']:
        report['warnings'].append('Input has exterior overhangs along this direction. Interface placement uses '
                                  'the outer ray envelope; thickness estimates can include voids. '
                                  'Inspect the exported layers before printing.')
    return [frame.world(p) for p in parts], report
