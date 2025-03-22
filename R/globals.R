# Declare global variables to avoid R CMD check warnings
# if (getRversion() >= "2.15.1") {
#   utils::globalVariables(c(
#     "ligand", "receptor", "sender", "receiver", "LR_pair", "SR_pair",
#     "cor", "p.scaled", "LRSR", "mean_score"
#   ))
# }
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "mean_score"
  ))
}
