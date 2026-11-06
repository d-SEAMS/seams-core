#!/usr/bin/env python3
"""Write the figure records that Snakemake keeps out of git.

Hydrate coordinates take oxygen from the GenIce dumps and place one
methane carbon at each cage centre. The nucleation host record and the
campaign manifest are the measurement records the statistics read.
"""
import base64
import io
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "paper" / "tex" / "figures" / "data"
XYZ = ROOT / "paper" / "tex" / "figures" / "xyzrender_xyz"

HOST = 'b41a6a061de80e941562fda8e4a18483395f72f646e799608a66689b88e1b02c  build/seams\nc3c75050bb7762cf0d3417e2fbab326422f23ea04a64b6f69eaf33334e4d1403  nucleation.lammpstrj\ncutoff identical yes\nseeded identical yes\nMEASURE2_EXIT=0\n'
MANIFEST = '{\n  "created": "2026-11-03T01:55:24Z",\n  "liquid_source": {\n    "name": "traj.lammpstrj",\n    "sha256": "d02023892336c43cf409da8ee4b7dee2b8f30f94ee9e76a1c1d1ae99e646ea2b",\n    "frame": 1\n  },\n  "confined_source": {\n    "name": "dump-240-square.lammpstrj",\n    "sha256": "ed256843b474ce06e724f6bb1caeba8e2a9f78b9aa6b5430048fe191e8790d77",\n    "doi": "10.6084/m9.figshare.11448768.v1",\n    "frame": 1,\n    "oxygen_type": 2\n  },\n  "frames": [\n    {\n      "system": "ic",\n      "scale": "1e3",\n      "molecules": 1000,\n      "box_A": [\n        31.75426480542942,\n        31.75426480542942,\n        31.75426480542942\n      ],\n      "npz": "ic-1e3.npz",\n      "trajectory": "ic-1e3.lammpstrj",\n      "npz_sha256": "051d248ca73f8548a99ec21ad871f246bd85e5470ff347232c8be0dd93fa8c16",\n      "trajectory_sha256": "79253e32eeaeb4a8e7b8166fe64e20690108ecbbdbc89ac3f9b392221ca8d4e4"\n    },\n    {\n      "system": "ic",\n      "scale": "1e4",\n      "molecules": 10648,\n      "box_A": [\n        69.85938257194472,\n        69.85938257194472,\n        69.85938257194472\n      ],\n      "npz": "ic-1e4.npz",\n      "trajectory": "ic-1e4.lammpstrj",\n      "npz_sha256": "0df8fde3790e9a3584fbc4c8e7a52bb4bbee7c0434072b7172ba1825c7a9c3e5",\n      "trajectory_sha256": "9659a909456d5c0dcf0637db637320158a9c31af3745d42b575308f494d2cfee"\n    },\n    {\n      "system": "ic",\n      "scale": "1e5",\n      "molecules": 97336,\n      "box_A": [\n        146.06961810497532,\n        146.06961810497532,\n        146.06961810497532\n      ],\n      "npz": "ic-1e5.npz",\n      "trajectory": "ic-1e5.lammpstrj",\n      "npz_sha256": "306e42a59f186f14f38bf68a1c4fe94782ccc13e94101bf1d5fa5f4d78782a94",\n      "trajectory_sha256": "f78f5c1ec05a3edffbc80b9ddc5f43b3eba29aa414f043ae58fb28bfe04d50b0"\n    },\n    {\n      "system": "ic",\n      "scale": "1e6",\n      "molecules": 1000000,\n      "box_A": [\n        317.54264805429415,\n        317.54264805429415,\n        317.54264805429415\n      ],\n      "npz": "ic-1e6.npz",\n      "trajectory": "ic-1e6.lammpstrj",\n      "npz_sha256": "8ceeadc37797dc6c02ba137ad42925f56865def3bb9883259fa1fee03fb694e7",\n      "trajectory_sha256": "c0b879d795405171b6686b65585c3850646a6ec02a606731be9c7b6fa85b0346"\n    },\n    {\n      "system": "ih",\n      "scale": "1e3",\n      "molecules": 1000,\n      "box_A": [\n        22.45365597551247,\n        38.890872965260115,\n        36.66666666666667\n      ],\n      "npz": "ih-1e3.npz",\n      "trajectory": "ih-1e3.lammpstrj",\n      "npz_sha256": "774905c1184c6e6a21cb491f0b52e54e207dae7e420a9fcefc7b0793e5d16c7b",\n      "trajectory_sha256": "a216733afe1584d5e621b06a3ae853b6563b554670000c3ac0f192eef5ea54b4"\n    },\n    {\n      "system": "ih",\n      "scale": "1e4",\n      "molecules": 10648,\n      "box_A": [\n        49.39804314612743,\n        85.55992052357225,\n        80.66666666666667\n      ],\n      "npz": "ih-1e4.npz",\n      "trajectory": "ih-1e4.lammpstrj",\n      "npz_sha256": "603a917ca3ccf3397eb8a402612b15e5b4bfde3158de96e76d275a2fd579ffaa",\n      "trajectory_sha256": "810cd14581fb9ed3618a26339d47c172233c63dd522d405d2f8270ca31bd3325"\n    },\n    {\n      "system": "ih",\n      "scale": "1e5",\n      "molecules": 97336,\n      "box_A": [\n        103.28681748735735,\n        178.89801564019652,\n        168.66666666666669\n      ],\n      "npz": "ih-1e5.npz",\n      "trajectory": "ih-1e5.lammpstrj",\n      "npz_sha256": "b530bd3e0db04f823d4d677abc15dfa471b121dea492ce70e7b51b061812bab6",\n      "trajectory_sha256": "6a7424f0104887820684a35134b655cf0c2b0ce9860fd3bb958d647a76555a51"\n    },\n    {\n      "system": "ih",\n      "scale": "1e6",\n      "molecules": 1000000,\n      "box_A": [\n        224.53655975512467,\n        388.90872965260115,\n        366.6666666666667\n      ],\n      "npz": "ih-1e6.npz",\n      "trajectory": "ih-1e6.lammpstrj",\n      "npz_sha256": "6568944a1d24d5cb8f266a8ea2199dd9af1e1541c2a7e667e4784458d1cdaa83",\n      "trajectory_sha256": "eb72e60d2fa7aba7aca81c7e69a0a99abe03e3b7fdffbadafe74a07bc252f228"\n    },\n    {\n      "system": "liquid",\n      "scale": "1e3",\n      "molecules": 4096,\n      "box_A": [\n        49.62270434822753,\n        49.62270434822753,\n        49.62270434822753\n      ],\n      "npz": "liquid-1e3.npz",\n      "trajectory": "liquid-1e3.lammpstrj",\n      "npz_sha256": "aaf92c93201b4b2376b7da37c173225ecf8dcbf865016c81d94fe05409ee6724",\n      "trajectory_sha256": "a41eec66c831ee277cbd20f17a21540e261fc91cbc691733980d7db933e75440"\n    },\n    {\n      "system": "liquid",\n      "scale": "1e4",\n      "molecules": 32768,\n      "box_A": [\n        99.24540869645506,\n        99.24540869645506,\n        99.24540869645506\n      ],\n      "npz": "liquid-1e4.npz",\n      "trajectory": "liquid-1e4.lammpstrj",\n      "npz_sha256": "6068bba795f403e70ead3bdbfb8e616c6efc50772f536b7e022bebc82b64d2a9",\n      "trajectory_sha256": "10ae1ef7611ed78283fecdc4fc6b86a3cb4acdd0629701a9589e9e92547795b4"\n    },\n    {\n      "system": "liquid",\n      "scale": "1e5",\n      "molecules": 110592,\n      "box_A": [\n        148.8681130446826,\n        148.8681130446826,\n        148.8681130446826\n      ],\n      "npz": "liquid-1e5.npz",\n      "trajectory": "liquid-1e5.lammpstrj",\n      "npz_sha256": "53420136683a40dd96c6116d2400de0fb944c6eec82e209d4a0b97e15b85f3c0",\n      "trajectory_sha256": "76bf0227dc761a4093ce162a9db627ec3d057c66f30da5b24f972145b067d82e"\n    },\n    {\n      "system": "liquid",\n      "scale": "1e6",\n      "molecules": 884736,\n      "box_A": [\n        297.7362260893652,\n        297.7362260893652,\n        297.7362260893652\n      ],\n      "npz": "liquid-1e6.npz",\n      "trajectory": "liquid-1e6.lammpstrj",\n      "npz_sha256": "5a5fffd3cd037fab23d302bf3267d33ba162c78da3e83e0dc88a9ac429eea633",\n      "trajectory_sha256": "dcf6019d3ef5d31effe16d3c8675123c73bc22b4ab06e152bb7e9e51472e2ce0"\n    },\n    {\n      "system": "confined",\n      "scale": "1e3",\n      "molecules": 1000,\n      "box_A": [\n        80.0,\n        80.0,\n        180.20265674635874\n      ],\n      "npz": "confined-1e3.npz",\n      "trajectory": "confined-1e3.lammpstrj",\n      "npz_sha256": "7678f2da78d7143bdad51fe98911c95f553c0cea3a474c136aa667b1e3ad642f",\n      "trajectory_sha256": "7cf6274f9b43ee1c91df2e1a45f2d78b3f6af6748df0cb23c8a9bcb07bcb6ac3"\n    },\n    {\n      "system": "confined",\n      "scale": "1e4",\n      "molecules": 10000,\n      "box_A": [\n        160.0,\n        200.0,\n        360.4053134927175\n      ],\n      "npz": "confined-1e4.npz",\n      "trajectory": "confined-1e4.lammpstrj",\n      "npz_sha256": "009bdbb69b8088f57179706cb58ca7c2f52b1994e28ef5877728db71139f101e",\n      "trajectory_sha256": "8258e4dd829c618aaa57b6b6bfd1dcc4caf7f451122a9927e064827ec9a0620b"\n    },\n    {\n      "system": "confined",\n      "scale": "1e5",\n      "molecules": 100000,\n      "box_A": [\n        400.0,\n        400.0,\n        720.810626985435\n      ],\n      "npz": "confined-1e5.npz",\n      "trajectory": "confined-1e5.lammpstrj",\n      "npz_sha256": "80235949a790db6265375ff78a3cb33082e82d95785f37cd6eeded16058bd139",\n      "trajectory_sha256": "4363a99ff33f9dd4779985e9d40966ecaf4a18b54614a7943fb68986abe62111"\n    },\n    {\n      "system": "confined",\n      "scale": "1e6",\n      "molecules": 1000000,\n      "box_A": [\n        800.0,\n        800.0,\n        1802.0265674635875\n      ],\n      "npz": "confined-1e6.npz",\n      "trajectory": "confined-1e6.lammpstrj",\n      "npz_sha256": "0cf4bc935acba67ac575ae4e6ac96eacb23c27bda6b059ff0649dfbf4fccabc1",\n      "trajectory_sha256": "003b1a285f8f2bbe5691298a5599d6ec2a33e08002ea1dde31858fb2cbdc1359"\n    }\n  ]\n}\n'
CAGES = {
    'si': [
        (6.0, 6.0, 6.0),
        (0.0, 0.0, 0.0),
        (0.0, 6.0, 9.0),
        (6.0, 9.0, 0.0),
        (6.0, 3.0, 0.0),
        (0.0, 6.0, 3.0),
        (9.0, 0.0, 6.0),
        (3.0, 0.0, 6.0),
    ],
    'sii': [
        (0.0, 0.0, 8.56),
        (12.84, 8.56, 12.84),
        (8.56, 12.84, 12.84),
        (8.56, 0.0, 0.0),
        (0.0, 12.84, 4.28),
        (0.0, 8.56, 0.0),
        (8.56, 8.56, 8.56),
        (8.56, 4.28, 4.28),
        (0.0, 4.28, 12.84),
        (12.84, 4.28, 0.0),
        (12.84, 0.0, 4.28),
        (4.28, 0.0, 12.84),
        (4.28, 4.28, 8.56),
        (4.28, 12.84, 0.0),
        (4.28, 8.56, 4.28),
        (12.84, 12.84, 8.56),
        (6.42, 6.42, 14.98),
        (14.98, 6.42, 6.42),
        (14.98, 14.98, 14.98),
        (6.42, 14.98, 6.42),
        (10.7, 10.7, 2.14),
        (10.7, 2.14, 10.7),
        (2.14, 2.14, 2.14),
        (2.14, 10.7, 10.7),
    ],
}


def lammps_frame(path):
    text = Path(path).read_text(encoding="utf-8").splitlines()
    index = 0
    while "NUMBER OF ATOMS" not in text[index]:
        index += 1
    count = int(text[index + 1])
    index += 2
    while not text[index].startswith("ITEM: BOX"):
        index += 1
    index += 1
    box = []
    while not text[index].startswith("ITEM:"):
        parts = text[index].split()
        box.append(float(parts[1]) - float(parts[0]))
        index += 1
    columns = text[index].split()[2:]
    index += 1
    points = []
    for row in text[index:index + count]:
        record = dict(zip(columns, row.split()))
        points.append((float(record["x"]), float(record["y"]), float(record["z"])))
    return points, box


def xyz_text(name, points, box, carbons):
    length = f"{box[0]:.8f} 0 0 0 {box[1]:.8f} 0 0 0 {box[2]:.8f}"
    lines = [
        str(len(points) + len(carbons)),
        f'{name} Lattice="{length}" Properties=species:S:1:pos:R:3',
    ]
    for x, y, z in points:
        lines.append(f" O {x:.8f} {y:.8f} {z:.8f}")
    for x, y, z in carbons:
        lines.append(f" C {x:.8f} {y:.8f} {z:.8f}")
    return "\n".join(lines) + "\n"


def write_hydrate():
    XYZ.mkdir(parents=True, exist_ok=True)
    specs = (
        ("si_hydrate", ROOT / "input/traj/genice_sI.lammpstrj", "si"),
        ("sii_hydrate", ROOT / "input/traj/genice_sII.lammpstrj", "sii"),
    )
    for name, source, key in specs:
        points, box = lammps_frame(source)
        (XYZ / f"{name}.xyz").write_text(xyz_text(name, points, box, CAGES[key]), encoding="utf-8")


def write_host():
    (DATA / "nucleation-host.txt").write_text(HOST, encoding="utf-8")


def write_manifest():
    path = DATA / "h2h" / "input-manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(MANIFEST, encoding="utf-8")


INPUTS_SHA = (
    "c3c75050bb7762cf0d3417e2fbab326422f23ea04a64b6f69eaf33334e4d1403  nucleation.lammpstrj\n"
    "ed256843b474ce06e724f6bb1caeba8e2a9f78b9aa6b5430048fe191e8790d77  dump-240-square.lammpstrj\n"
    "d02023892336c43cf409da8ee4b7dee2b8f30f94ee9e76a1c1d1ae99e646ea2b  liquid.lammpstrj\n"
)


def write_records():
    from recorded_inputs import PAYLOAD

    raw = base64.b64decode("".join(PAYLOAD.split()))
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        archive.extractall(ROOT, filter="data")
    (DATA / ".records-ready").write_text("ok\n", encoding="utf-8")
    sha_path = ROOT / "paper" / "tex" / "repro" / "inputs.sha256"
    sha_path.write_text(INPUTS_SHA, encoding="utf-8")


def main():
    command = sys.argv[1] if len(sys.argv) > 1 else "all"
    if command in ("host", "all"):
        write_host()
    if command in ("manifest", "all"):
        write_manifest()
    if command in ("hydrate", "all"):
        write_hydrate()
    if command in ("records", "all"):
        write_records()


if __name__ == "__main__":
    main()
