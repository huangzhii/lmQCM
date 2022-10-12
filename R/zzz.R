# first and last lib functions

.onAttach = function(libname, pkgname)
{
  ourVer = try( gsub("[^0-9_.-]", "", utils::packageVersion("lmQCM"), fixed = FALSE) );

  if (inherits(ourVer, "try-error")) ourVer = "";
  packageStartupMessage(paste("  Package lmQCM", ourVer, "loaded.\n"))
  packageStartupMessage("==========================================================================\n");
  packageStartupMessage(paste0(
         " If you benefit from this package, please cite:\n",
         " \n",
         "Huang Z, Han Z, Wang T, Shao W, Xiang S, Salama P, Rizkalla M, Huang K, Zhang J.\n",
         "TSUNAMI: Translational Bioinformatics Tool Suite For Network Analysis And Mining.\n",
         "Genomics, Proteomics & Bioinformatics, Volume 19, Issue 6, December 2021, Pages 1023-1031\n",
         "https://doi.org/10.1016/j.gpb.2019.05.006."))
  packageStartupMessage("==========================================================================\n\n");
}

