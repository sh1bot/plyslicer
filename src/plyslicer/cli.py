import argparse
import json
import sys
from pathlib import Path

import numpy as np
import trimesh

from .geometry import GeometryError, load_model, split_model


def vector(value):
    try:
        result = [float(x) for x in value.split(',')]
        if len(result) != 3:
            raise ValueError
        return result
    except ValueError:
        raise argparse.ArgumentTypeError('Expected X,Y,Z') from None


def main(argv=None):
    parser = argparse.ArgumentParser(description='Split a thin STL/3MF solid into curved, straight-assembling laminations.')
    parser.add_argument('input', type=Path)
    parser.add_argument('-o', '--output', type=Path, required=True, help='New output directory (never overwrites)')
    parser.add_argument('-n', '--layers', type=int, default=2)
    parser.add_argument('--pitch', type=float, default=.8, help='Interface sampling spacing in mm (default: .8)')
    parser.add_argument('--direction', type=vector, help='Override automatic insertion direction: X,Y,Z in world coordinates')
    parser.add_argument('--min-thickness', type=float, default=.4, help='Minimum 5th-percentile axial layer thickness in mm')
    parser.add_argument('--pins', type=int, default=0, help='Number of removable alignment-rod bores, 0–8')
    parser.add_argument('--pin-diameter', type=float, default=2., help='Alignment rod diameter in mm')
    parser.add_argument('--clearance', type=float, default=.15, help='Radial rod clearance in mm')
    parser.add_argument('--print-oriented', action='store_true', help='Also export alternate build directions in a print/ subdirectory; supports may be needed')
    args = parser.parse_args(argv)
    try:
        if args.output.exists():
            raise GeometryError('Output directory already exists; choose a new directory.')
        if args.input.suffix.lower() not in ('.stl', '.3mf'):
            raise GeometryError('Input must be STL (coordinates assumed mm) or 3MF.')
        print('Loading and validating mesh...', file=sys.stderr, flush=True)
        mesh = load_model(args.input)
        print('Selecting direction, sampling curved interfaces, and cutting layers...', file=sys.stderr, flush=True)
        parts, report = split_model(mesh, args.layers, args.pitch, args.direction,
                                    args.min_thickness, args.pins, args.pin_diameter, args.clearance)
        args.output.mkdir(parents=True)
        for i, part in enumerate(parts, 1):
            part.export(args.output / f'layer-{i:02d}.stl')
        if args.print_oriented:
            directory = args.output / 'print'
            directory.mkdir()
            local = np.eye(4)
            local[:3, :3] = np.array(report['frame_basis_columns']).T
            local[:3, 3] = -local[:3, :3] @ np.array(report['frame_origin_world'])
            for i, part in enumerate(parts):
                rotation = trimesh.transformations.rotation_matrix(-np.pi/2 if i % 2 else 0, [0, 1, 0])
                transform = rotation @ local
                printable = part.copy()
                printable.apply_transform(transform)
                shift = np.eye(4)
                shift[:2, 3] = -printable.bounds.mean(axis=0)[:2]
                shift[2, 3] = -printable.bounds[0, 2]
                printable.apply_transform(shift)
                printable.export(directory / f'layer-{i+1:02d}.stl')
                report['parts'][i]['assembly_to_print_transform'] = (shift @ transform).tolist()
        report['input'] = str(args.input)
        (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        print(f'Exported {len(parts)} closed layers to {args.output}')
        print(f'Bulk layer thickness (5th percentile): {report["bulk_layer_thickness_p05_mm"]:.3f} mm')
        print('Read report.json for assembly direction, alignment, and print limitations.')
        return 0
    except (GeometryError, OSError, ValueError) as exc:
        print(f'plyslicer: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
