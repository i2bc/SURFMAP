
import setuptools


def get_license(): 
    with open('LICENSE') as f:
        license = f.read()

    return license


packages = [
    'bin',
    'surfmap',
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
                'surfmap=bin.surfmap:main',
                'extract_interface=bin.extract_interface:main',
                '_surfmap_tool=surfmap.tools.SurfmapTools:main'
                ]
            }
    )