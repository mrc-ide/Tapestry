import os
from typing import List, Dict

# TODO
# - Sanity checks and useful warnings
# - Could change to using Path

class TapestryCOIOutputDir:
	"""
	Interface for `Tapestry` output directory for
	a particular COI
	"""
	def __init__(self, coi_dir: str) -> None:
		self.coi_dir = coi_dir

		# MCMC diagnostics
		self.mcmc_diagnostics_csv = f"{coi_dir}/mcmc.diagnostics.csv"
		self.mcmc_parameters_csv = f"{coi_dir}/mcmc.parameters.csv"

		# Fits
		self.props_csv = f"{coi_dir}/fit.proportions.csv"
		self.ibd_pairwise_csv = f"{coi_dir}/fit.ibd.pairwise.csv"
		self.ibd_path_csv = f"{coi_dir}/fit.ibd.path.csv"
		self.ibd_segs_bed = f"{coi_dir}/fit.ibd.segments.bed"
		self.ibd_stats_csv = f"{coi_dir}/fit.ibd.stats.csv"
	

class TapestrySampleOutputDir:
	"""
	Interface for outputs produced by `Tapestry` for
	a single sample

	"""

	def __init__(self, tapestry_dir: str) -> None:
		"""
		Initialise interface for all input files

		"""
		
		self.tapestry_dir = tapestry_dir

		# COI Information
		self._coi_dirs, self._cois = self._get_coi_info()
		
		# COI dictionary -- inteface for interacting with COI specific files
		self.coi = self._build_coi_dict()

		# Comparison CSVs
		self.compare_heuristic_csv = f"{self.tapestry_dir}/compare.heuristic.csv"
		self.compare_evidence_csv = f"{self.tapestry_dir}/compare.evidence.csv"

	def _get_coi_info(self) -> (List[str], List[int]):
		""" 
		Get all COI directories, 'K[0-9]{1}'
		
		"""

		coi_dirs = [
			f"{self.tapestry_dir}/{d}"
			for d in os.listdir(self.tapestry_dir)
			if d.startswith("K") and os.path.isdir(f"{self.tapestry_dir}/{d}")
		]
		cois = [int(coi_dir[-1]) for coi_dir in coi_dirs]
		
		return (coi_dirs, cois)
	
	def _build_coi_dict(self) -> Dict[int, TapestryCOIOutputDir]:
		"""
		Build a dictionary of COI outputs

		"""
		assert self._coi_dirs is not None
		assert self._cois is not None
		return {
			coi: TapestryCOIOutputDir(coi_dir)
			for coi, coi_dir in zip(self._cois, self._coi_dirs)
		}
	
