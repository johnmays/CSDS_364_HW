# Notes for the feedback survey for A3:

### A3a:

2b. `gensignal` is confusing. It says it should take arguments `(t, g, tau, T)` but then it specifies an example in which `fs` is given as a parameter. Additionally, with respect to the example, it is impossible for the function `plot_sampled_function(g, ...)` to plot `gensignal()` because `gensignal()` must be initialized with both `t` and another function `g`. When I pass `plot_sampled_function(g, ...)` a function such as `sin()`, I am essentially passing it the function name alone, and then calling it with the `t` that is created within, however that would be impossible to do here because I need to call it on whatever function will constitute the signal. So now it is unclear to me whether or not I need to write a new plotting function. The instructions on how to write the function and what it must do are fine I suppose, but how to "show" it? Very confusing to me.

For 3b and 2b, I spent a lot of time questioning whether or not these functions are going to be passed continuous time values or not. Perhaps it could be made more clear which will suit the rest of the assignment so I don't have to rewrite later.
