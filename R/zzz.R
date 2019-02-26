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
         " Zhang, Jie & Huang, Kun (2014) \n",
         " Normalized ImQCM: An Algorithm for Detecting Weak Quasi-Cliques\n",
         " in Weighted Graph with Applications in Gene Co-Expression Module\n",
         " Discovery in Cancers. Cancer informatics, 13, CIN-S14021.\n"))
  packageStartupMessage("==========================================================================\n\n");
}

