from pathlib import Path
import os
import shutil
import subprocess
import tempfile


class MosaicRunner(object):
    def __init__(
        self, simulation_name=None, root=".", cpl_dir="/usr/src", user_name="root"
    ):
        if simulation_name is None:
            simulation_name = next(tempfile._get_candidate_names()).replace("x", "_")
        else:
            if len(simulation_name) > 15:
                raise AttributeError(
                    f"simulation name: {simulation_name} is too long (>15)"
                )
        self.simulation_name = simulation_name
        self.root = root
        self.cpl_dir = cpl_dir
        self.user_name = user_name  # used to send data at end to /home/user_name

        # TODO These should probably be but in child classes
        self.AtmMdl = "LMD96x95"
        self.OceMdl = "PALEORCA2"

        self.run_def = DefaultRunConfig()
        self.run_def.c_basins = ".".join([self.OceMdl, simulation_name, "nc"])
        self.run_def.cname = ".".join([self.OceMdl, simulation_name])

    @property
    def domsk_dir(self):
        return os.path.join(
            self.root, "ORCA2_BLD", "DOMSK", self.OceMdl + "." + self.simulation_name
        )

    @property
    def mosaic_dir(self):
        return os.path.join(
            self.root,
            "ORCA2_BLD",
            self.OceMdl + "." + self.simulation_name + "x" + self.AtmMdl,
        )

    def create_directory_structure(self):
        # Create DOMSK directory
        Path(self.domsk_dir).mkdir(parents=True, exist_ok=True)
        Path(self.mosaic_dir).mkdir(parents=True, exist_ok=True)

    def cleanup(self):
        shutil.rmtree(self.domsk_dir)
        shutil.rmtree(self.mosaic_dir)

    def setup_domsk(self, bathy_file, coordinates_file, subbasins_file):
        # Create run def
        run_def_filepath = self.domsk_dir + "/run.def"
        self.run_def.to_file(run_def_filepath)

        # Copy files to directory with correctname.
        shutil.copy2(bathy_file, self.domsk_dir + "/bathy_meter.nc")
        shutil.copy2(coordinates_file, self.domsk_dir + "/coordinates.nc")
        shutil.copy2(
            subbasins_file, self.domsk_dir + f"/subbasins_{self.run_def.cname}.nc"
        )

    def run_domsk(self):
        cmd = os.path.join(self.cpl_dir, "DOMSK", "bin", "domsk.exe")
        proc = subprocess.Popen(cmd, cwd=self.domsk_dir)
        # Wait for process to finish
        proc.wait()
        # rename output files

    def run_file_generation(self):
        # Setup Architecture
        cmd = os.path.join(self.cpl_dir, "util", "do_link.sh")
        proc = subprocess.Popen(cmd, cwd=self.mosaic_dir).wait()

        # Add run.def
        mozaic_dir = os.path.join(self.mosaic_dir, "MOZAIC")
        run_def_filepath = os.path.join(mozaic_dir, "run.def")
        self.run_def.to_file(run_def_filepath)

        # Run Grids.sh
        print("Running Grids", flush=True)
        subprocess.Popen("./grids.sh", cwd=mozaic_dir, stdout=subprocess.PIPE).wait()

        print("Running mosaic", flush=True)
        proc = subprocess.Popen("./mosaic.exe", cwd=mozaic_dir, stdout=subprocess.PIPE)
        (output, err) = proc.communicate()
        proc.wait()

        print("Updating rundef jma2o and jmo2a")
        output = output.decode("utf-8").split("\n")
        for line in output:
            if "oce -> atm ADRESSE1 WEIGHTS1 Neighbors :" in line:
                _continue = True
                jmo2a = line.rstrip().lstrip().split()[-1]
                while _continue:
                    try:
                        int(jmo2a[-1])
                        # Stop as soon as we hit a digit
                        _continue = False
                    except ValueError:
                        # Remove last character
                        jmo2a = jmo2a[:-1]
            if "atm -> oce ADRESSE2 WEIGHTS2 Number of neighbors :" in line:
                _continue = True
                jma2o = line.rstrip().lstrip().split()[-1]
                while _continue:
                    try:
                        int(jma2o[-1])
                        # Stop as soon as we hit a digit
                        _continue = False
                    except ValueError:
                        # Remove last character
                        jma2o = jma2o[:-1]

        # Update run def values
        self.run_def.jmo2a = jmo2a
        self.run_def.jma2o = jma2o
        self.run_def.to_file(run_def_filepath)

        # Run Grids.sh -r
        print("Running Grids", flush=True)
        subprocess.Popen(
            ["./grids.sh", "-r"], cwd=mozaic_dir, stdout=subprocess.PIPE
        ).wait()

        # Running cotes
        print("Running cotes", flush=True)
        proc = subprocess.Popen("./cotes.exe", cwd=mozaic_dir, stdout=subprocess.PIPE)
        (output, err) = proc.communicate()
        proc.wait()

        print("Updating rundef jma2or")
        output = output.decode("utf-8").split("\n")
        for line in output:
            if "atm -> oce ADRESSE3 WEIGHTS3 Number of neighbors :" in line:
                _continue = True
                jma2or = line.rstrip().lstrip().split()[-1]
                while _continue:
                    try:
                        int(jma2or[-1])
                        # Stop as soon as we hit a digit
                        _continue = False
                    except ValueError:
                        # Remove last character
                        jma2or = jma2or[:-1]

        self.run_def.jma2or = jma2or
        self.run_def.to_file(run_def_filepath)

        # Run allwei
        print("Running allwei")
        subprocess.Popen(
            ["./allwei.sh"], cwd=self.mosaic_dir, stdout=subprocess.PIPE
        ).wait()

        # Copy out needed files
        subprocess.Popen(
            ["./envoie.sh", "-i", "-5A2", "-L", "39", "-D", "-u", self.user_name],
            cwd=self.mosaic_dir,
            stdout=subprocess.PIPE,
        ).wait()

    def run(self, bathy_file, coords_file, subbasins_file):
        self.create_directory_structure()
        self.setup_domsk(bathy_file, coords_file, subbasins_file)
        self.run_domsk()
        self.run_file_generation()


class RunConfig(object):
    def to_file(self, filepath="run.def"):
        with open(filepath, "w") as f:
            for key, value in vars(self).items():
                if isinstance(value, (tuple, list)):
                    value = ", ".join([str(v) for v in value])
                f.write(str(key) + " = " + str(value) + "\n")


class DefaultRunConfig(RunConfig):
    def __init__(self):
        # IPSL en mode debug ?
        self.l_ipsldbg = "n"
        #
        self.c_suffix = "none"
        #
        # Debug allocs ?
        self.l_d_alloc = "n"
        #
        # Dryrun (test I/O) ?
        self.l_dryrun = "n"
        self.lev_dry = 0
        #
        # Cas LGM ou non
        self.clgm = ".FALSE."
        # config pour domsk.exe
        self.l_paleorca2 = "y"
        self.clcoo = "coordinates.nc"
        self.cname = "paleorca2"
        self.cl_met = "bathy_meter.nc"
        self.cr = "rp2"
        self.l_bat_cdf = "y"
        self.group = "PALEORCA2"
        self.production = "NEMO - PALEORCA2"
        self.jpi = 182
        self.jpj = 149
        self.nperio = 6
        self.jpk = 31
        self.l_meter = "y"
        self.l_wri_3D = "y"
        #
        ## Modele ocean
        # Dimensions
        self.jpoi = 182
        self.jpoj = 149
        self.l_recalc_o = "n"
        # Type de periodicite
        self.noperio = 6
        # Nombre de points pour decrire la maille
        self.jpoe = 9
        # Nom du modele
        self.comod = "orc"
        # Type de la grille
        self.cotyp = "orca2.3"
        # Ordre des noms
        self.locerev = "n"
        #
        ## Modele atmosphere
        # Dimensions
        self.jpai = 96
        self.jpaj = 95
        self.l_recalc_a = "n"
        # Nombre de points pour decrire la maille
        self.jpae = 9
        # Type de periodicite
        self.naperio = -1
        self.la_pole = "y"
        self.la_nortop = "n"
        # Nom du modele
        self.camod = "lmd"
        # Type de la grille
        self.catyp = "lmdz"

        # Type of normalization atm->oce : 0: none, 1: intensive, 2: extensive
        self.norma2o = 2
        self.normo2a = 1

        # Max number of weights (dimension of arrays)
        self.jpa2o = 10000
        self.jpo2a = 100
        #
        # Actual number of weights
        # o -> a
        self.jmo2a = 15
        # a -> o fluxes
        self.jma2o = 18
        # a -> o runoff
        self.jma2or = 9
        # a -> o calving
        self.jma2oi = 50
        #
        # Parametres calcul de poids de runoff
        # Traitement du runoff des rivieres avec les embouchures exactes
        self.lriv = "n"
        # Traitement specifique des points cotiers
        self.lcoast = "y"
        # Calcul pour run-off integre sur la maille atm
        self.lint_atm = "y"
        # Calcul pour run-off integre sur la maille oce
        self.lint_oce = "n"
        # Extension de 1 point a l''interieur, vers le point ocean le plus proche
        self.lnear = "n"
        # Extension de 1 point a l''interieur, vers le point atm voisin
        self.lnei = "y"
        # Tout les points atmospheres mouilles les plus proches.
        self.ltotal = "n"
        # Tout les points atmospheres mouilles les plus proches. Avec distance maximum
        self.ltotal_dist = "y"
        self.dist_max = 400.0e3
        # Tout les points atmospheres mouilles les plus proches. Avec distance maximum et etalement sur voisins proches ocean
        self.ltotal_dist_2 = "y"
        self.dist_max_voisin = 500.0e3

        # What ocean mask is written in Mosaic MCT file for runoff : noperio/perio
        self.cotes_omsk = "noperio"

        # What atmosphere mask is written in Mosaic MCT file for runoff : int/ext/full
        self.cotes_amsk = "full"

        # Le masque ocean du run-off est sans les points doubles
        self.l_runoff_msk_noperio = "y"  ## Non actif
        # Le masque atmosphere du run-off est extensif (tous les points avec iun fraction de terre sont comptes terre)
        self.l_runoff_msk_ext = "y"  ## Non actif

        # Parametre de calcul des poids pour le calvin

        # Fichier du masque de bassin
        # En fait egal a paleorca2.nc qui est sorti dans DOMSK lorsque l_wri_3D = y
        self.c_basins = "PALEORCA2.20MaF.nc"  ##eORCA1.2.nc
        # Nom des variables
        self.cl_atl = "mask_atl"
        self.cl_pac = "mask_pac"
        self.cl_nomed = "mask_nomed"
        self.cl_noclo = "mask_noclose"

        # Nombre de bandes de latitude
        self.jp_calv = 3
        # Limites des bandes
        self.ylimits = -90.0, -55.0, 40.0, 90.0
        # Suppression du Pacific nord
        self.l_calving_nopac = "n"
        # Suppression de l'Atlantique nord
        self.l_calving_noatl = "n"
        # Suppression des Mediterrannee (inclus Mer Rouge, Golfe Persique)
        self.l_calving_nomed = "y"

        #
        ##
        # Use limit stack to reduce memory in NetCDF (for big cases only)
        self.limit_stack = ".FALSE."
        # Type of OASIS (2.2 for all version after)
        self.c_oasis = 2.2
        # Use if OASIS grid are available in NetCDF format
        self.l_grid_cdf = ".TRUE."
        # Format NetCDF (syntaxe Fliocom) 'REPLACE', 'REPHDF'
        self.c_FlioMode = "REPHDF"
        #
        self.l_wei_i4 = "y"
        self.l_wei_i8 = "n"
        self.l_wei_oasis3 = "y"
        self.l_wei_oasis_mct = "y"


if __name__ == "__main__":
    r = MosaicRunner("test")
    r.create_directory_structure()
    r.setup_domsk(
        "/tmp/data/bathyORCA2.RupelianTotalV1.nc",
        "/tmp/data/coordinates_paleorca2_40Ma.nc",
    )
    r.run_domsk()
    r.run_file_generation()
