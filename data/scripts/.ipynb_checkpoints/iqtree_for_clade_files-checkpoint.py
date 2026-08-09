import subprocess
from pathlib import Path

folders = [
    "../Aves_Amphibia_trees/P_Y_clade_parts_2"
]

for folder in folders:
    folder = Path(folder)

    print(f"Обработка папки: {folder}")

    for file in folder.glob("*.fasta"):
        print(f"Запуск IQ-TREE для: {file}")

        cmd = [
            "iqtree",
            "-s", str(file),
            "-st", "DNA",
            "-m", "MFP",
            "-bb", "1000",
            "-nt", "AUTO"
        ]

        subprocess.run(cmd, check=True)

print("Все выравнивания обработаны")