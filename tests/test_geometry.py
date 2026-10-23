import json

import numpy as np
import pytest
import trimesh

from plyslicer.cli import main
from plyslicer.geometry import GeometryError, load_model, make_frame, projection_check, slab, split_model


def curved_shell():
    x, y = np.meshgrid(np.linspace(-10, 10, 21), np.linspace(-20, 20, 41))
    bottom = .015*x*x + .002*y*y
    return slab(np.column_stack((x.ravel(), y.ravel())), bottom, bottom + 3)


@pytest.mark.parametrize('layers', [2, 3, 5])
def test_flat_closed_volume_and_thickness(layers):
    mesh = trimesh.creation.box([20, 40, 5])
    parts, report = split_model(mesh, layers=layers, pitch=1, direction=[0, 0, 1])
    assert len(parts) == layers
    assert all(p.is_watertight for p in parts)
    assert [p.volume for p in parts] == pytest.approx([mesh.volume/layers]*layers, rel=1e-5)
    assert report['bulk_layer_thickness_p05_mm'] == pytest.approx(5/layers)


def test_curved_interfaces_and_straight_assembly():
    mesh = curved_shell()
    parts, report = split_model(mesh, layers=3, pitch=.7, direction=[0, 0, 1])
    assert sum(p.volume for p in parts) == pytest.approx(mesh.volume, rel=1e-5)
    assert all(p.is_watertight and len(p.split()) == 1 for p in parts)
    # A curved interface must contain significant height variation.
    assert np.ptp(parts[0].vertices[:, 2]) > 2
    for shift in [.1, 1, 10]:
        upper = parts[1].copy()
        upper.apply_translation([0, 0, shift])
        overlap = trimesh.boolean.intersection([parts[0], upper], engine='manifold')
        assert overlap.is_empty or overlap.volume < 1e-5


def test_exterior_undercut_is_not_an_interface_undercut():
    a = trimesh.creation.box([20, 20, 2])
    b = a.copy()
    b.apply_translation([0, 0, 6])
    bridge = trimesh.creation.box([2, 20, 8])
    bridge.apply_translation([9, 0, 3])
    mesh = trimesh.boolean.union([a, b, bridge], engine='manifold')
    check = projection_check(make_frame(mesh, [0, 0, 1]).local(mesh))
    assert not check['passes']
    parts, _ = split_model(mesh, direction=[0, 0, 1], pitch=1)
    upper = parts[1].copy()
    upper.apply_translation([0, 0, .1])
    overlap = trimesh.boolean.intersection([parts[0], upper], engine='manifold')
    assert overlap.is_empty or overlap.volume < 1e-5


def test_alignment_bores_remove_material_and_keep_closed_parts():
    mesh = trimesh.creation.box([30, 60, 4])
    parts, report = split_model(mesh, pitch=1, pins=2, direction=[0, 0, 1])
    assert len(report['alignment']['bore_centers_local_xy_mm']) == 2
    assert 0 < sum(p.volume for p in parts) < mesh.volume
    assert all(p.is_watertight for p in parts)


def test_cli_exports_and_refuses_overwrite(tmp_path):
    path = tmp_path / 'box.stl'
    trimesh.creation.box([10, 20, 4]).export(path)
    out = tmp_path / 'parts'
    args = [str(path), '-o', str(out), '--print-oriented']
    assert main(args) == 0
    report = json.loads((out / 'report.json').read_text())
    assert report['layers'] == 2
    for file in out.glob('*.stl'):
        assert trimesh.load_mesh(file).is_watertight
    axes = []
    for part in report['parts']:
        transform = np.array(part['assembly_to_print_transform'])
        axes.append(transform[:3, :3].T @ [0, 0, 1])
        printable = trimesh.load_mesh(out / 'print' / f'layer-{part["layer"]:02d}.stl')
        assert printable.bounds[0, 2] == pytest.approx(0, abs=1e-5)
        assert printable.is_watertight
    assert axes[0] @ axes[1] == pytest.approx(0, abs=1e-8)
    assert main(args) == 2


def test_open_mesh_rejected(tmp_path):
    mesh = trimesh.creation.box()
    mesh.update_faces(np.arange(len(mesh.faces)-1))
    path = tmp_path / 'open.stl'
    mesh.export(path)
    assert main([str(path), '-o', str(tmp_path / 'out')]) == 2


def test_3mf_units_and_build_transform(tmp_path):
    import io
    import zipfile
    import xml.etree.ElementTree as ET

    scene = trimesh.Scene()
    transform = np.eye(4)
    transform[:3, 3] = [3, 4, 5]
    scene.add_geometry(trimesh.creation.box([1, 2, .2]), transform=transform)
    data = scene.export(file_type='3mf')
    path = tmp_path / 'inches.3mf'
    with zipfile.ZipFile(io.BytesIO(data)) as source, zipfile.ZipFile(path, 'w') as dest:
        for name in source.namelist():
            content = source.read(name)
            if name.endswith('.model'):
                root = ET.fromstring(content)
                root.set('unit', 'inch')
                content = ET.tostring(root)
            dest.writestr(name, content)
    mesh = load_model(path)
    assert mesh.extents == pytest.approx(np.array([1, 2, .2])*25.4)
    assert mesh.bounds.mean(axis=0) == pytest.approx(np.array([3, 4, 5])*25.4)


@pytest.mark.parametrize('kwargs', [{'layers': 1}, {'pitch': 0}, {'direction': [0, 0, 0]},
                                    {'pins': -1}, {'clearance': -1}, {'pitch': float('nan')}])
def test_bad_parameters(kwargs):
    with pytest.raises(GeometryError):
        split_model(trimesh.creation.box([10, 20, 4]), **kwargs)
