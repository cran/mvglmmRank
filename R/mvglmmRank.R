mvglmmRank <-
function (game.data, method = "PB0", first.order = FALSE, home.field=TRUE, max.iter.EM = 1000, 
    tol1 = 1e-04, tol2 = 1e-04, tolFE = 0, tol.n = 1e-07, verbose = TRUE,OT.flag=FALSE) 
{
    if (class("method") != "character") {
        cat("*Error: method must be a character string (using quotation marks)")
        flush.console()
        return(0)
    }
    if (first.order) 
        control <- list(iter.EM = max.iter.EM, tol1 = tol1, tol2 = tol2, 
            tolFE = tolFE, verbose = verbose,OT.flag=OT.flag)
    if (!first.order) 
        control <- list(iter.EM = max.iter.EM, tol1 = tol1, tol2 = tol2, 
            tolFE = tolFE, verbose = verbose,OT.flag=OT.flag)
    control.n <- list(iter.EM = max.iter.EM, tol1 = tol.n, verbose = verbose,OT.flag=OT.flag)
    Z_mat <- game.data
    Z_mat$Score.For <- Z_mat$home.score
    Z_mat$Score.Against <- Z_mat$away.score
    Z_mat$home_win <- as.numeric(Z_mat$Score.For > Z_mat$Score.Against)
    Z_mat$home <- as.character(Z_mat$home)
    Z_mat$away <- as.character(Z_mat$away)
    if (is.null(Z_mat$neutral.site)) 
        Z_mat$neutral.site <- 0
    if (method == "B") {
        res <- binary_cre(Z_mat = Z_mat, first.order = first.order, home.field=home.field, 
            control = control)
    }
    else if (method == "P0") {
        res <- poisson_cre(Z_mat = Z_mat, first.order = first.order,home.field=home.field, 
            control = control, game.effect = FALSE)
    }
    else if (method == "P1") {
        res <- poisson_cre(Z_mat = Z_mat, first.order = first.order,home.field=home.field, 
            control = control, game.effect = TRUE)
    }
    else if (method == "N") {
        res <- normal_cre(Z_mat = Z_mat, first.order = first.order,home.field=home.field, 
            control = control.n)
    }
    else if (method == "NB") {
        res <- NB_cre(Z_mat = Z_mat, first.order = first.order, home.field=home.field,
            control = control)
    }
    else if (method == "PB0") {
        res <- PB_cre(Z_mat = Z_mat, first.order = first.order, home.field=home.field, 
            control = control, game.effect = FALSE)
    }
    else if (method == "PB1") {
        res <- PB_cre(Z_mat = Z_mat, first.order = first.order,  home.field=home.field,
            control = control, game.effect = TRUE)
    }
    else {
        cat("Error in specification of method. This field is case sensitive.\n")
        return(0)
    }
    class(res) <- "mvglmmRank"
    return(res)
}
