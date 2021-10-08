#' @title Main variational EM algorithm for SLDS
#' @description Generate approximate posterior distribution of continous and discrete latent
#'     state variables and point estimates for model parameterslibrary("benchmarkme")
#' @param y observations
#' @param qs matrix of size 2 by TT with each row being the probability vector for S_t
#' @param k2 observational variance when st = 2
#' @param is.obs is h_mat obseved value or likelihood value
#'
#'
#' @return list of variational parameters of q(X) and q(S), ELBO, and q(X) and q(S)
#' @export
#'


variational_em_slds <- function(y, qs = matrix(0.5, 2, TT),
                                k2 = 30, is.obs = TRUE) {

    ################################
    # Arguments

    ################################

    TT <- length(y)
    # ################################
    # Initialization of parameters
    # ################################

    ## Assume model parameters are KNOWN at this time,
    ## and focus on E-step to see if VI part works
    phi <- 0.9
    v <- 1
    w <- 0.5
    xi <- 0.99
    xi_mat <- matrix(c(xi, 1 - xi, 1 - xi, xi), 2, 2)

    ## {F, G, V, W, m0, C0} under two different states.
    ## In this example, only V is state-dependent
    FF_lst <- list(1, 1)
    GG_lst <- list(phi, phi)
    V_lst <- list(v, k2 * v)
    W_lst <- list(w, w)
    m0_lst <- list(0, 0)
    C0_lst <- list(1, 1)
    V_inv_lst <- lapply(V_lst, function(x) 1 / x)

    d <- 1 ## dimension of continous latent states
    ELBO_old <- 1

    # filtering storage
    a_mat <- matrix(0, nrow = d, ncol = TT)
    R_arry <- array(0, dim = c(d, d, TT))
    m_mat <- matrix(0, nrow = d, ncol = TT)
    C_arry <- array(0, dim = c(d, d, TT))

    # smoothing storage
    ms_mat <- matrix(0, nrow = d, ncol = TT)
    Cs_arry <- array(0, dim = c(d, d, TT))


    # objects used often
    V_inv_FF_lst <- Map('%*%', V_inv_lst, FF_lst)
    W_inv_lst <- lapply(W_lst, function(x) {chol2inv(chol(x))})
    FF_trans <- lapply(FF_lst, function(x) t(x))
    FF_V_inv_FF_lst <- Map('%*%', FF_trans, V_inv_FF_lst)
    GG_trans <- lapply(GG_lst, function(x) t(x))
    W_inv_GG_lst <- Map('%*%', W_inv_lst, GG_lst)
    GG_W_inv_GG_lst <- Map('%*%', GG_trans, W_inv_GG_lst)
    C0_inv_lst <- lapply(C0_lst, function(x) {chol2inv(chol(x))})
    C0_inv_m0_lst <- Map('%*%', C0_inv_lst, m0_lst)


    ## Do not run while loop this time becasue the ELBO is NOT increasing and the
    ## algorithm does not converge

    # while(eps_em > epsilon_em) {

        # if (count == 100) {
        #     print("Does not converge")
        #     break
        # }

        # ################################
        # ******** E-step ************ #
        # ################################
    esp_elbo <- 1

    ## iterate 100 times
    ELBO_vec <- rep(0, 100)
        # ==============================
        # Variational approximation
        # ==============================
    for (mm in 1:100) {
        cat("mm = ", mm, "\r")

        Vhat_inv_lst <- vector("list", TT)
        What_inv_lst <- vector("list", TT)
        FFhat_lst <- vector("list", TT)
        GGhat_lst <- vector("list", TT)

        for (t in TT:1) {
            qs_t <- qs[, t]

            # -------------
            ### update Vhat
            # -------------
            Vhat_inv <- update_Vhat_inv(V_inv_lst, qs_t)
            Vhat_inv_lst[[t]] <- Vhat_inv           ## save V_hat_inv

            # -------------
            ### update FFhat
            # -------------
            FFhat <- update_FFhat(V_inv_FF_lst, Vhat_inv, qs_t)
            FFhat_lst[[t]] <- FFhat           ## save FFhat

            # -------------
            ### update What_inv
            # -------------
            if (t >= 2) {
                if (t == TT) {
                    What_inv <- update_What_inv(W_inv_lst = W_inv_lst,
                                                FF_V_inv_FF_lst = FF_V_inv_FF_lst,
                                                Vhat_inv = Vhat_inv,
                                                FFhat = FFhat,
                                                qs_t = qs_t, t = t, TT = TT)
                } else {
                    What_inv <- update_What_inv(W_inv_lst = W_inv_lst,
                                                FF_V_inv_FF_lst = FF_V_inv_FF_lst,
                                                Vhat_inv = Vhat_inv,
                                                FFhat = FFhat,
                                                GG_W_inv_GG_lst = GG_W_inv_GG_lst,
                                                What_inv = What_inv_lst[[t + 1]],
                                                GGhat = GGhat_lst[[t + 1]],
                                                qs_t = qs_t, t = t, TT = TT)
                }
                What_inv_lst[[t]] <- What_inv
            }


            # -------------
            ### update Ghat_inv
            # -------------
            if (t >= 2) {
                GGhat <- update_GGhat(What_inv, W_inv_GG_lst, qs_t)
                GGhat_lst[[t]] <- GGhat
            }

            # -------------
            ### update C0hat
            # -------------
            if (t == 1) {
                C0hat_inv <- update_C0hat_inv(C0_inv_lst, FF_V_inv_FF_lst,
                                              Vhat_inv, FFhat,
                                              GG_W_inv_GG_lst = GG_W_inv_GG_lst,
                                              What_inv = What_inv_lst[[t + 1]],
                                              GGhat = GGhat_lst[[t + 1]],
                                              qs_t)
            }

            # -------------
            ### update m0hat
            # -------------
            if(t == 1) {
                m0hat <- update_m0hat(C0hat_inv, C0_inv_m0_lst, qs_t)
            }
        }

        # ====================================================
        # Update q(x | lambda_x) by Kalman smoothing
        # ====================================================
        Vhat_lst <- lapply(Vhat_inv_lst, function(x) 1 / x)
        What_lst <- lapply(What_inv_lst, function(x) 1 / x)
        C0hat <- chol2inv(chol(C0hat_inv))

        Vhat_mat <- matrix(unlist(Vhat_lst), TT, 1)
        mod <- dlm::dlm(FF = FFhat, V = 1, JV = 1,
                   GG = GGhat, W = 1,
                   m0 = m0hat, C0 = C0hat,
                   X = Vhat_mat)
        smooth_res <- dlm::dlmSmooth(dlm::dlmFilter(y, mod))
        Cs_vec <- unlist(dlm::dlmSvd2var(smooth_res$U.S, smooth_res$D.S))[2:(TT + 1)]
        Cs <- array(Cs_vec, dim = c(1, 1, TT))
        # Cs <- array(dlmSvd2var(smooth_res$U.S, smooth_res$D.S), dim = c(1, 1, 200))
        ms <- matrix(smooth_res$s[2:(TT + 1)], 1, TT)

        ## My own code
        # # ---------
        # # filtering
        # # ---------
        # for (t in 2:TT) {
        #     if (t == 2) {
        #         at <- GGhat_lst[[t]] %*% m0hat
        #         Rt <- GGhat_lst[[t]] %*% C0hat %*% t(GGhat_lst[[t]]) + What_lst[[t]]
        #     } else {
        #         at <- GGhat_lst[[t]] %*% m_mat[, t - 1]
        #         Rt <- GGhat_lst[[t]] %*% C_arry[, , t - 1] %*% t(GGhat_lst[[t]]) + What_lst[[t]]
        #     }
        #     et <- y[t] - FFhat_lst[[t]] %*% at
        #     RFt <- Rt %*% t(FFhat_lst[[t]])
        #     Qt <- FFhat_lst[[t]] %*% RFt + Vhat_lst[[t]]
        #     At <- RFt %*% chol2inv(chol(Qt))
        #     mt <- at + At %*% et
        #     Ct <- Rt - At %*% Qt %*% t(At)
        #     a_mat[, t] <- at
        #     R_arry[, , t] <- Rt
        #     m_mat[, t] <- mt
        #     C_arry[, , t] <- Ct
        # }
        # # ---------
        # # smoothing
        # # ---------
        # ms_mat[, TT] <- m_mat[, TT]
        # Cs_arry[, , TT] <- C_arry[, , TT]
        #
        # for (t in (TT - 1):1) {
        #     Bt <- C_arry[, , t] %*%
        #         t(GGhat_lst[[t + 1]]) %*% chol2inv(chol(R_arry[, , t + 1]))
        #     ms <- m_mat[, t] + Bt %*% (ms_mat[, t + 1] - a_mat[, t + 1])
        #     Cs <- C_arry[, , t] +
        #         Bt %*% (Cs_arry[, , t + 1] - R_arry[, , t + 1]) %*% t(Bt)
        #     ms_mat[, t] <- ms
        #     Cs_arry[, , t] <- Cs
        # }
        # ms <- ms_mat
        # Cs <- Cs_arry
        # # if(mm == 5) {
        # #     ms_5 <- ms_mat
        # #     Cs_5 <- Cs_arry
        # # }
        # Cs[, , 1] <- 1

        # ====================================================
        # Compute parameters of discrete latent states S (2 states)
        # ====================================================
        log_h_mat <- matrix(0, nrow = 2, ncol = TT)

        ### only for this constant case
        FF <- FF_lst[[1]]
        GG <- GG_lst[[1]]
        W <- W_lst[[1]]
        m0 <- m0_lst[[1]]
        C0 <- C0_lst[[1]]

        ### t = 1
        for (k in 1:2) {
            V <- V_lst[[k]]
            V_inv <- 1 / V
            V_inv_FF <- V_inv %*% FF
            FF_V_inv_FF <- emulator::quad.form(V_inv, FF)

            B <- V_inv %*% y[1] ^ 2  - 2 * V_inv_FF %*% m0hat * y[1] +
                FF_V_inv_FF %*% (C0hat + m0hat ^ 2)

            Tr <- sum(diag(B))
            log_det_V <- log(V)
            log_h_mat[k, 1] <- -1/2 * (Tr + log_det_V) - 1/2 - 1/2 * log(C0hat)
        }


        ### t = 2:TT
        for (t in 2:TT) {

            Exx <- sufficient_stats(ms, Cs, tm1 = t, t = t)
            ## replace GGhat_lst[[t]] by phi for AR(1) case
            Ex1x <- sufficient_stats(ms, Cs, tm1 = t - 1, t = t,
                                     GG = GGhat_lst[[t]], m0, C0)
            Ex1x1 <- sufficient_stats(ms, Cs, tm1 = t - 1, t = t - 1)

            # W_inv <- chol2inv(W)
            W_inv <- 1 / W
            W_inv_GG <- W_inv %*% GG
            GG_W_inv_GG <- emulator::quad.form(W_inv, GG)
            A <- W_inv %*% Exx - 2 * W_inv_GG %*% Ex1x + GG_W_inv_GG %*% Ex1x1

            for (k in 1:2) {
                V <- V_lst[[k]]
                # V_inv <- chol2inv(V)
                V_inv <- 1 / V
                V_inv_FF <- V_inv %*% FF
                FF_V_inv_FF <- emulator::quad.form(V_inv, FF)

                B <- V_inv %*% y[t] ^ 2  - 2 * V_inv_FF %*% ms[, t] * y[t] +
                    FF_V_inv_FF %*% Exx

                Tr <- sum(diag(A + B))
                # log_det_W <- determinant(W, logarithm = TRUE)$modulus
                # log_det_V <- determinant(V, logarithm = TRUE)$modulus
                log_det_W <- log(W)
                log_det_V <- log(V)

                log_h_mat[k, t] <- -1/2 * (Tr + log_det_W + log_det_V)
            }
        }


        # ====================================================
        # Update q(s | lambda_s) by forward-backward algorithm for HMM
        # ====================================================

        h_mat <- exp(log_h_mat)
        h_mat_normal <- apply(h_mat, 2, function(x) x/sum(x))
        if (is.obs) {
            ### h_mat as observed value
            fit_hmm <- HiddenMarkov::Estep(x = h_mat_normal[1, ],
                                           Pi = xi_mat, delta = c(0.5, 0.5),
                                           distn = "beta",
                                           pm = list(shape1 = c(1, 8),
                                                     shape2 = c(8, 1)))
            qs <- t(fit_hmm$u)  ## (TT x 2)

        } else {            ### h_mat as likelihood value
            qs <- state_prob_hmm(K = 2, n = TT, h_mat = h_mat,
                                 init_dist = c(0.5, 0.5), tran_prob = xi_mat)

            qs_sum <- apply(qs, 2, sum)
            for (i in 1:TT) {
                qs[, i] <- qs[, i] / qs_sum[i]
            }
        }

        # ====================================================
        # Compute ELBO
        # ====================================================
        ELBO_y <- 0
        ELBO_x <- 0
        ELBO_s <- 0  ## constant
        ELBO_qx <- 0
        ELBO_qs <- 0

        ELBO_x1 <- -1/2 * ((Cs[, , 1] + ms[, 1] ^ 2) - 2 * ms[, 1] * m0 + m0 ^ 2) / C0

        for (t in 1:TT) {
            avg_V <- wt_mean(V_lst, qs[, t])
            avg_logV <- wt_mean(lapply(V_lst, log), qs[, t])
            V_inv <- lapply(V_lst, function(x) 1 / x)
            avg_V_inv <- wt_mean(V_inv, qs[, t])
            Exx <- sufficient_stats(ms, Cs, tm1 = t, t = t)

            Ex1x <- sufficient_stats(ms, Cs, tm1 = t - 1, t = t, GG = phi,
                                     m0 = m0, C0 = C0)
            Ex1x1 <- sufficient_stats(ms, Cs, tm1 = t - 1, t = t - 1)

            ELBO_y <- ELBO_y + avg_logV +
                avg_V_inv * (y[t] ^ 2 - 2 * y[t] * ms[, t] + Exx)

            if (t >= 2) {
                ELBO_x <- ELBO_x + Exx - 2 * phi * Ex1x + phi ^ 2 * Ex1x1
            }

            ELBO_qs <- ELBO_qs + stats::weighted.mean(log(qs[, t]), qs[, t])
        }

        ELBO_y <- -1/2 * ELBO_y
        ELBO_x <- -TT/2 * log(W) - 1/2 * ELBO_x / W
        ELBO_qx <- -1/2 * apply(log(Cs), c(1, 2), sum) - TT / 2

        (ELBO <- ELBO_y + ELBO_x + ELBO_x1 + ELBO_s - ELBO_qx - ELBO_qs)
        ELBO_vec[mm] <- ELBO
        (esp_elbo <- ELBO - ELBO_old)

        ELBO_old <- ELBO
    }

        # ################################
        # ******** M-step ************ #
        # ################################
        ### Assume model parameters are known and fixed at this time.


    # }
    return(list(qs = qs,
                log_h_mat = log_h_mat,
                ms = ms,
                Cs = Cs,
                Vhat_inv_lst = Vhat_inv_lst,
                What_inv_lst = What_inv_lst,
                FFhat_lst = FFhat_lst,
                GGhat_lst = GGhat_lst,
                C0hat = C0hat,
                m0hat = m0hat,
                ELBO_vec = ELBO_vec))
}
