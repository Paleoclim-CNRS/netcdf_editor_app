import os
import glob
import shutil
import subprocess
import tempfile

# from coords_gen import coordinates_generator


class MosaixRunner(object):

    def __init__(self, bathy_file, coords_file):
        self.mosaix_dir = os.path.join('/', 'src', 'MOSAIX')
        self.coordinates_mask = None
        self.bathy_file = bathy_file # Filename
        self.coords_file = coords_file # Filename
        self.output_dir = next(tempfile._get_candidate_names())
        if not os.path.exists(self.output_dir):
            os.mkdir(os.path.join(self.mosaix_dir, self.output_dir))

    def create_coordinates_mask(self):
        # coordiantes_generator(self.bathy_file, self.coords_file)
        # Create Coordinates mask
        # set self.coordiantes_mask with the current filename (not necessarily in right location)
        self.coordinates_mask = "/root/Scratch/IGCM/OCE/NEMO/ORCA2.3/ORCA2.3_coordinates_mask.nc"

    def run_mosaix(self):
        if self.coordinates_mask is None:
            raise AttributeError("Coordinates mask is not set please run self.create_coordinates_mask()")
        # Copy coordinates_xios.nc to coordiantes_mask.nc into Mosaix
        try:
            shutil.copy2(self.coordinates_mask, os.path.join(os.environ['HOME'], 'Scratch', 'IGCM', 'OCE', 'NEMO', 'ORCA2.3'))
        except shutil.SameFileError:
            pass
        proc = subprocess.Popen("./CreateWeightsMask.bash", cwd=self.mosaix_dir)
        # Wait for process to finish
        proc.wait()

        files_to_copy = [
            os.path.join(self.mosaix_dir, 'README_ORCA2.3xLMD9695_MOSAIX_v1.txt'),
            *glob.glob(os.path.join(self.mosaix_dir, 'areas_ORCA2.3*')),
            *glob.glob(os.path.join(self.mosaix_dir, 'grids_ORCA2.3*')),
            *glob.glob(os.path.join(self.mosaix_dir, 'masks_ORCA2.3*')),
            *glob.glob(os.path.join(self.mosaix_dir, 'LMD9695_grid_maskFrom_ORCA2.3*.nc')),
            *glob.glob(os.path.join(self.mosaix_dir, 'dia_*ORCA2.3*.nc')),
            *glob.glob(os.path.join(self.mosaix_dir, 'rmp_*ORCA2.3*.nc')),
        ] 

        for _file in files_to_copy:
            print(f"Moving file {_file} to {os.path.join(self.mosaix_dir, self.output_dir, os.path.basename(_file))}", flush=True)
            shutil.move(
                _file, 
                os.path.join(
                    self.mosaix_dir, self.output_dir, os.path.basename(_file)
                    )
                )


    def run(self):
        self.create_coordinates_mask()
        self.run_mosaix()

    def cleanup(self):
        shutil.rmtree(os.path.join(self.mosaix_dir, self.output_dir))