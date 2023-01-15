
import setuptools


def get_license(): 
    with open('LICENSE') as f:
        license = f.read()

    return license


packages = [
    'surfmap',
    'surfmap.bin',
    'surfmap.tools',
    'surfmap.lib',
    ]

setuptools.setup(
    name='surfmap',
    python_requires='>=3.7',
    license=get_license(),
    packages=packages,
    include_package_data=True,
    install_requires=[
        "pdb2pqr==3.1.0",
        "requests==2.22.0",
        "numpy==1.24.1",
        "freesasa==2.1.0",
        ],
        entry_points={
            'console_scripts': [
                'surfmap=surfmap.bin.surfmap:main',
                '_surfmap_tool=surfmap.tools.SurfmapTools:main',
                'extract_interface=surfmap.bin.extract_interface:main',
                'write_pdb_bs=surfmap.bin.write_pdb_interface:main',
                ]
            }
    )