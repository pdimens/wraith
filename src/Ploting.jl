using VegaLite

function plotJaccard(mat::Array{Int64,2}, filename::String)
    p = mat |>
    @vlplot(
        :rect,
        width = 300, height = 200,
        x = {:IMDB_Rating, bin = {maxbins = 60}},
        y = {:Rotten_Tomatoes_Rating, bin = {maxbins = 40}},
        color = "count()",
        config = {
            range = {
                heatmap = {
                    scheme = "greenblue"
                }
            },
            view = {
                stroke = "transparent"
            }
        }
    )
    if filename
        save(filename, p)
    end
    return p
end