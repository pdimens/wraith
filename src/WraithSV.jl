module WraithSV

include("Jaccard.jl")

export jaccardIdent, jaccardScores, DiagonalStats, diagonalStats, diagonalSD,
    gompertz, gompertzJacobianRow, fitGompertz, predictionBands, findOutliers,
    detectOutliers

end
