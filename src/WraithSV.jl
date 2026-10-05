module WraithSV

include("Jaccard.jl")
include("PackedTriangle.jl")

export jaccardIdent, jaccardScores, DiagonalStats, diagonalStats, diagonalSD,
    gompertz, gompertzJacobianRow, fitGompertz, predictionBands, findOutliers,
    detectOutliers, PackedTriangle, packedlength, packTriangle, unpack,
    writePacked, readPacked, openMapped

end
