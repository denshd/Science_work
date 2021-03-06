using PlotlyJS
# Build the rectangles as a heatmap
# specify the edges of the heatmap squares
phi = (1 + sqrt(5) )/2. # golden ratio
xe = [0, 1, 1+(1/(phi^4)), 1+(1/(phi^3)), phi]
ye = [0, 1/(phi^3), 1/phi^3+1/phi^4, 1/(phi^2), 1]

z = [13 3 3 5; 13 2 1 5; 13 10 11 12; 13 8 8 8]'

trace = heatmap(
    x=sort(xe),
    y=sort(ye),
    z=z,
    type="heatmap",
    colorscale="Viridis"
)
plot(trace)