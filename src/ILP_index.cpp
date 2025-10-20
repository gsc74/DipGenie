#include "ILP_index.h"

// Constructor
ILP_index::ILP_index(gfa_t *g) : Solver(g) {}

void printQuadraticConstraints(GRBModel& model) {
    GRBQConstr* qconstraints = model.getQConstrs();
    int numQConstraints = model.get(GRB_IntAttr_NumQConstrs);

    for (int i = 0; i < numQConstraints; ++i) {
        std::string qconstrName = qconstraints[i].get(GRB_StringAttr_QCName);
        GRBQuadExpr qconstrExpr = model.getQCRow(qconstraints[i]);
        double qrhs = qconstraints[i].get(GRB_DoubleAttr_QCRHS);
        char qsense = qconstraints[i].get(GRB_CharAttr_QCSense);

        std::cout << "Quadratic Constraint " << qconstrName << ": ";

        // Print linear terms
        for (int j = 0; j < qconstrExpr.getLinExpr().size(); ++j) {
            GRBVar var = qconstrExpr.getLinExpr().getVar(j);
            double coeff = qconstrExpr.getLinExpr().getCoeff(j);
            std::string varName = var.get(GRB_StringAttr_VarName);

            std::cout << coeff << " * " << varName;
            if (j < qconstrExpr.getLinExpr().size() - 1 || qconstrExpr.size() > 0) {
                std::cout << " + ";
            }
        }

        // Print quadratic terms
        for (int j = 0; j < qconstrExpr.size(); ++j) {
            GRBVar var1 = qconstrExpr.getVar1(j);
            GRBVar var2 = qconstrExpr.getVar2(j);
            double qcoeff = qconstrExpr.getCoeff(j);
            std::string varName1 = var1.get(GRB_StringAttr_VarName);
            std::string varName2 = var2.get(GRB_StringAttr_VarName);

            std::cout << qcoeff << " * " << varName1 << " * " << varName2;
            if (j < qconstrExpr.size() - 1) {
                std::cout << " + ";
            }
        }

        if (qsense == GRB_EQUAL) {
            std::cout << " == ";
        } else if (qsense == GRB_LESS_EQUAL) {
            std::cout << " <= ";
        } else if (qsense == GRB_GREATER_EQUAL) {
            std::cout << " >= ";
        }

        std::cout << qrhs << std::endl;
    }

    delete[] qconstraints;
}

void printConstraints(GRBModel& model) {
    GRBConstr* constraints = model.getConstrs();
    int numConstraints = model.get(GRB_IntAttr_NumConstrs);

    for (int i = 0; i < numConstraints; ++i) {
        std::string constrName = constraints[i].get(GRB_StringAttr_ConstrName);
        GRBLinExpr constrExpr = model.getRow(constraints[i]);
        double rhs = constraints[i].get(GRB_DoubleAttr_RHS);
        char sense = constraints[i].get(GRB_CharAttr_Sense);

        std::cout << "Constraint " << constrName << ": ";

        for (int j = 0; j < constrExpr.size(); ++j) {
            GRBVar var = constrExpr.getVar(j);
            double coeff = constrExpr.getCoeff(j);
            std::string varName = var.get(GRB_StringAttr_VarName);

            std::cout << coeff << " * " << varName;
            if (j < constrExpr.size() - 1) {
                std::cout << " + ";
            }
        }

        if (sense == GRB_EQUAL) {
            std::cout << " == ";
        } else if (sense == GRB_LESS_EQUAL) {
            std::cout << " <= ";
        } else if (sense == GRB_GREATER_EQUAL) {
            std::cout << " >= ";
        }

        std::cout << rhs << std::endl;
    }

    delete[] constraints;
}

void printNonZeroVariables(GRBModel& model) {
    int numVars = model.get(GRB_IntAttr_NumVars);
    GRBVar* vars = model.getVars();

    std::cout << "Non-zero variables:\n";
    for (int i = 0; i < numVars; ++i) {
        double val = vars[i].get(GRB_DoubleAttr_X);
        if (val != 0.0) {
            std::string varName = vars[i].get(GRB_StringAttr_VarName);
            std::cout << varName << " = " << val << std::endl;
        }
    }
    delete[] vars;
}

void printObjectiveFunction(GRBModel& model) {
    GRBQuadExpr objExpr = model.getObjective();
    int sense = model.get(GRB_IntAttr_ModelSense);

    std::cout << "Objective Function: ";

    // Print linear terms
    bool firstTerm = true; // To handle the '+' sign placement
    for (int i = 0; i < objExpr.getLinExpr().size(); ++i) {
        GRBVar var = objExpr.getLinExpr().getVar(i);
        double coeff = objExpr.getLinExpr().getCoeff(i);
        std::string varName = var.get(GRB_StringAttr_VarName);

        if (!firstTerm) {
            std::cout << " + ";
        }
        std::cout << coeff << " * " << varName;
        firstTerm = false;
    }

    // Print quadratic terms
    for (int i = 0; i < objExpr.size(); ++i) {
        GRBVar var1 = objExpr.getVar1(i);
        GRBVar var2 = objExpr.getVar2(i);
        double qcoeff = objExpr.getCoeff(i);
        std::string varName1 = var1.get(GRB_StringAttr_VarName);
        std::string varName2 = var2.get(GRB_StringAttr_VarName);

        if (!firstTerm) {
            std::cout << " + ";
        }
        std::cout << qcoeff << " * " << varName1 << " * " << varName2;
        firstTerm = false;
    }

    // Print constant term if it exists
    double constant = model.get(GRB_DoubleAttr_ObjCon);
    if (constant != 0) {
        if (!firstTerm) {
            std::cout << " + ";
        }
        std::cout << constant;
    }

    if (sense == GRB_MINIMIZE) {
        std::cout << " (Minimize)" << std::endl;
    } else {
        std::cout << " (Maximize)" << std::endl;
    }
}

void ILP_index::solve(std::vector<std::pair<std::string, std::string>> &ip_reads)
{



    // For backtracking the haplotype
    std::map<int, std::string> haplotypes; // h -> haplotype string
    // Write an ILP with Gurobi
    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        // env.set(GRB_IntParam_Threads, 0);
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Set parameters to speed up the model
        model.set("PreSparsify", "1"); // Sparsify the model
        model.set("Heuristics", "0.50"); // Spent 50% time on heuristics
        model.set("NodefileStart", "0.5"); // 0.5 GB nodefile start
        model.set("Presolve", "2"); // Aggressive presolve to reduce the model size
        model.set("Method", "3"); // Concurrent method
        model.set("Threads", std::to_string(num_threads));   // let Gurobi decide

        // Initialize the objective function and expressions
        GRBLinExpr obj;
        GRBLinExpr vtx_expr;
        GRBLinExpr z_expr; // For the objective function

        // Data structures
        std::map<int, std::map<std::string, GRBVar>> vars;          // h -> var_name -> GRBVar
        std::map<int, std::map<int32_t, GRBVar>> Zvars_h_alpha;           // h -> i -> GRBVar
        std::map<int, std::map<int32_t, GRBVar>> Zvars_h_beta;           // h -> i -> GRBVar
        std::map<int32_t, GRBVar> Zvars_alpha;                            // i -> GRBVar
        std::map<int32_t, GRBVar> Zvars_beta;                            // i -> GRBVar

        int32_t c_1 = recombination_penalty; // INF no recombination

        bool is_ilp = true; // ILP
        if (is_qclp) is_ilp = false; // QCLP
        int32_t count_kmer_matches = 0;

        if (is_ilp)
        {
            fprintf(stderr, "[M::%s::%.3f*%.2f] ILP model started\n",
                    __func__, realtime() - mg_realtime0, cputime()/(realtime()-mg_realtime0));

            // ── Clear all per-haplotype maps exactly once ─────────────────
            for (int h = 1; h <= ploidy; ++h) {
                vars[h].clear();
                Zvars_h_alpha[h].clear();
                Zvars_h_beta[h].clear();
            }

            // ── Homozygous K-mer constraints ─────────────────────────────
            for (int h = 1; h <= ploidy; ++h) {
                for (int32_t i = 0; i < count_sp_r; i++) {
                    GRBLinExpr z_expr_h;
                    int32_t temp = 0;

                    for (int32_t j = 0; j < num_walks; j++) {
                        if (paths[j].empty()) continue;
                        for (int32_t k = 0; k < Anchor_hits_homo[i][j].size(); k++) {
                            GRBLinExpr kmer_expr;
                            std::string extra_var   = "alpha_" + std::to_string(i) + "_" +
                                                    std::to_string(j) + "_" +
                                                    std::to_string(k);
                            std::string extra_var_h = extra_var + "_" + std::to_string(h);
                            GRBVar kmer_expr_var    = model.addVar(
                                0.0, 1.0, 0.0, GRB_BINARY, extra_var_h
                            );

                            if (Anchor_hits_homo[i][j][k].size() - 1 == 0)
                                continue;

                            for (int32_t l = 1; l < Anchor_hits_homo[i][j][k].size(); l++) {
                                int32_t u = Anchor_hits_homo[i][j][k][l-1];
                                int32_t v = Anchor_hits_homo[i][j][k][l];
                                std::string var_name   = std::to_string(u) + "_" +
                                                        std::to_string(j) + "_" +
                                                        std::to_string(v) + "_" +
                                                        std::to_string(j);
                                std::string var_name_h = var_name + "_" + 
                                                        std::to_string(h);

                                if (vars[h].find(var_name) == vars[h].end()) {
                                    vars[h][var_name] = model.addVar(
                                        0.0, 1.0, 0.0,
                                        is_mixed ? GRB_CONTINUOUS : GRB_BINARY,
                                        var_name_h
                                    );
                                }
                                kmer_expr += vars[h][var_name];
                            }

                            int32_t weight = Anchor_hits_homo[i][j][k].size() - 1;
                            model.addConstr(
                                kmer_expr >= weight * kmer_expr_var,
                                "Kmer_constraints_alpha_" + std::to_string(i) + "_" +
                                std::to_string(j) + "_" + std::to_string(k) + "_" +
                                std::to_string(h)
                            );

                            z_expr_h += kmer_expr_var;
                            temp++;
                        }
                    }

                    if (temp != 0) {
                        std::string z_var   = "alpha_" + std::to_string(i);
                        std::string z_var_h = z_var + "_" + std::to_string(h);
                        GRBVar z_var_r      = model.addVar(
                            0.0, 1.0, 0.0, GRB_BINARY, z_var_h
                        );
                        Zvars_h_alpha[h][i] = z_var_r;
                        model.addConstr(
                            z_expr_h == z_var_r,
                            "alpha_constraint_" + std::to_string(i) + "_" +
                            std::to_string(h)
                        );
                        if (h == 1) count_kmer_matches++;
                    }
                }
            }

            // ── Heterozygous K-mer constraints ────────────────────────────
            for (int h = 1; h <= ploidy; ++h) {
                for (int32_t i = 0; i < count_sp_r; i++) {
                    GRBLinExpr z_expr_h;
                    int32_t temp = 0;

                    for (int32_t j = 0; j < num_walks; j++) {
                        if (paths[j].empty()) continue;
                        for (int32_t k = 0; k < Anchor_hits_hetero[i][j].size(); k++) {
                            GRBLinExpr kmer_expr;
                            std::string extra_var   = "beta_" + std::to_string(i) + "_" +
                                                    std::to_string(j) + "_" +
                                                    std::to_string(k);
                            std::string extra_var_h = extra_var + "_" + std::to_string(h);
                            GRBVar kmer_expr_var    = model.addVar(
                                0.0, 1.0, 0.0, GRB_BINARY, extra_var_h
                            );

                            if (Anchor_hits_hetero[i][j][k].size() - 1 == 0)
                                continue;

                            for (int32_t l = 1; l < Anchor_hits_hetero[i][j][k].size(); l++) {
                                int32_t u = Anchor_hits_hetero[i][j][k][l-1];
                                int32_t v = Anchor_hits_hetero[i][j][k][l];
                                std::string var_name   = std::to_string(u) + "_" +
                                                        std::to_string(j) + "_" +
                                                        std::to_string(v) + "_" +
                                                        std::to_string(j);
                                std::string var_name_h = var_name + "_" +
                                                        std::to_string(h);

                                if (vars[h].find(var_name) == vars[h].end()) {
                                    vars[h][var_name] = model.addVar(
                                        0.0, 1.0, 0.0,
                                        is_mixed ? GRB_CONTINUOUS : GRB_BINARY,
                                        var_name_h
                                    );
                                }
                                kmer_expr += vars[h][var_name];
                            }

                            int32_t weight = Anchor_hits_hetero[i][j][k].size() - 1;
                            model.addConstr(
                                kmer_expr >= weight * kmer_expr_var,
                                "Kmer_constraints_beta_" + std::to_string(i) + "_" +
                                std::to_string(j) + "_" + std::to_string(k) + "_" +
                                std::to_string(h)
                            );

                            z_expr_h += kmer_expr_var;
                            temp++;
                        }
                    }

                    if (temp != 0) {
                        std::string z_var   = "beta_" + std::to_string(i);
                        std::string z_var_h = z_var + "_" + std::to_string(h);
                        GRBVar z_var_r      = model.addVar(
                            0.0, 1.0, 0.0, GRB_BINARY, z_var_h
                        );
                        Zvars_h_beta[h][i] = z_var_r;
                        model.addConstr(
                            z_expr_h == z_var_r,
                            "beta_constraint_" + std::to_string(i) + "_" +
                            std::to_string(h)
                        );
                        if (h == 1) count_kmer_matches++;
                    }
                }
            }
        }
        else
        {
            fprintf(stderr, "[M::%s::%.3f*%.2f] QP model started\n",
                    __func__, realtime() - mg_realtime0, cputime()/(realtime()-mg_realtime0));

            // ── Clear all per-haplotype maps exactly once ─────────────────
            for (int h = 1; h <= ploidy; ++h) {
                vars[h].clear();
                Zvars_h_alpha[h].clear();
                Zvars_h_beta[h].clear();
            }

            // ── Homozygous Q-constraints ────────────────────────────────
            for (int h = 1; h <= ploidy; ++h) {
                for (int32_t i = 0; i < count_sp_r; i++) {
                    GRBQuadExpr kmer_expr;
                    GRBLinExpr   z_expr_h;
                    int32_t      temp = 0;

                    for (int32_t j = 0; j < num_walks; j++) {
                        if (paths[j].empty()) continue;
                        for (int32_t k = 0; k < Anchor_hits_homo[i][j].size(); k++) {
                            std::string extra_var   = "alpha_" + std::to_string(i) + "_" +
                                                    std::to_string(j) + "_" +
                                                    std::to_string(k);
                            std::string extra_var_h = extra_var + "_" + std::to_string(h);
                            GRBVar   z_var_q        = model.addVar(
                                0.0, 1.0, 0.0, GRB_BINARY, extra_var_h
                            );

                            if (Anchor_hits_homo[i][j][k].size() - 1 == 0)
                                continue;

                            int32_t weight = Anchor_hits_homo[i][j][k].size() - 1;
                            for (int32_t l = 1; l < Anchor_hits_homo[i][j][k].size(); l++) {
                                int32_t u = Anchor_hits_homo[i][j][k][l-1];
                                int32_t v = Anchor_hits_homo[i][j][k][l];
                                std::string var_name   = std::to_string(u) + "_" +
                                                        std::to_string(j) + "_" +
                                                        std::to_string(v) + "_" +
                                                        std::to_string(j);
                                std::string var_name_h = var_name + "_" +
                                                        std::to_string(h);

                                if (vars[h].find(var_name) == vars[h].end()) {
                                    vars[h][var_name] = model.addVar(
                                        0.0, 1.0, 0.0,
                                        is_mixed ? GRB_CONTINUOUS : GRB_BINARY,
                                        var_name_h
                                    );
                                }
                                kmer_expr += vars[h][var_name] * z_var_q;
                            }
                            kmer_expr += (1 - weight) * z_var_q;

                            z_expr_h += z_var_q;
                            temp++;
                        }
                    }

                    if (temp != 0) {
                        std::string constraint_name = "Kmer_constraints_alpha_" +
                                                    std::to_string(i) + "_" +
                                                    std::to_string(h);
                        std::string z_var   = "alpha_" + std::to_string(i);
                        std::string z_var_h = z_var + "_" + std::to_string(h);
                        GRBVar   z_var_r    = model.addVar(
                            0.0, 1.0, 0.0, GRB_BINARY, z_var_h
                        );
                        Zvars_h_alpha[h][i] = z_var_r;
                        model.addQConstr(
                            kmer_expr == z_var_r,
                            constraint_name
                        );
                        model.addConstr(
                            z_expr_h == z_var_r,
                            "alpha_constraint_" + std::to_string(i) + "_" +
                            std::to_string(h)
                        );
                        if (h == 1) count_kmer_matches++;
                    }
                }
            }

            // ── Heterozygous Q-constraints ───────────────────────────────
            for (int h = 1; h <= ploidy; ++h) {
                for (int32_t i = 0; i < count_sp_r; i++) {
                    GRBQuadExpr kmer_expr;
                    GRBLinExpr   z_expr_h;
                    int32_t      temp = 0;

                    for (int32_t j = 0; j < num_walks; j++) {
                        if (paths[j].empty()) continue;
                        for (int32_t k = 0; k < Anchor_hits_hetero[i][j].size(); k++) {
                            std::string extra_var   = "beta_" + std::to_string(i) + "_" +
                                                    std::to_string(j) + "_" +
                                                    std::to_string(k);
                            std::string extra_var_h = extra_var + "_" + std::to_string(h);
                            GRBVar   z_var_q        = model.addVar(
                                0.0, 1.0, 0.0, GRB_BINARY, extra_var_h
                            );

                            if (Anchor_hits_hetero[i][j][k].size() - 1 == 0)
                                continue;

                            int32_t weight = Anchor_hits_hetero[i][j][k].size() - 1;
                            for (int32_t l = 1; l < Anchor_hits_hetero[i][j][k].size(); l++) {
                                int32_t u = Anchor_hits_hetero[i][j][k][l-1];
                                int32_t v = Anchor_hits_hetero[i][j][k][l];
                                std::string var_name   = std::to_string(u) + "_" +
                                                        std::to_string(j) + "_" +
                                                        std::to_string(v) + "_" +
                                                        std::to_string(j);
                                std::string var_name_h = var_name + "_" +
                                                        std::to_string(h);

                                if (vars[h].find(var_name) == vars[h].end()) {
                                    vars[h][var_name] = model.addVar(
                                        0.0, 1.0, 0.0,
                                        is_mixed ? GRB_CONTINUOUS : GRB_BINARY,
                                        var_name_h
                                    );
                                }
                                kmer_expr += vars[h][var_name] * z_var_q;
                            }
                            kmer_expr += (1 - weight) * z_var_q;

                            z_expr_h += z_var_q;
                            temp++;
                        }
                    }

                    if (temp != 0) {
                        std::string constraint_name = "Kmer_constraints_beta_" +
                                                    std::to_string(i) + "_" +
                                                    std::to_string(h);
                        std::string z_var   = "beta_" + std::to_string(i);
                        std::string z_var_h = z_var + "_" + std::to_string(h);
                        GRBVar   z_var_r    = model.addVar(
                            0.0, 1.0, 0.0, GRB_BINARY, z_var_h
                        );
                        Zvars_h_beta[h][i] = z_var_r;
                        model.addQConstr(
                            kmer_expr == z_var_r,
                            constraint_name
                        );
                        model.addConstr(
                            z_expr_h == z_var_r,
                            "beta_constraint_" + std::to_string(i) + "_" +
                            std::to_string(h)
                        );
                        if (h == 1) count_kmer_matches++;
                    }
                }
            }
        }

        // Combine Zvars_h[1][i] and Zvars_h[2][i] to create Zvars[i]
        std::set<int32_t> z_indices_alpha;
        std::set<int32_t> z_indices_beta;

        for (int h = 1; h <= ploidy; ++h) {
            for (const auto& z_pair : Zvars_h_alpha[h]) {
                z_indices_alpha.insert(z_pair.first);
            }
            for (const auto& z_pair : Zvars_h_beta[h]) {
                z_indices_beta.insert(z_pair.first);
            }
        }

        for (int32_t i : z_indices_alpha) {
            std::string z_var_name_alpha = "alpha_" + std::to_string(i);
            GRBVar z_var_alpha = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, z_var_name_alpha);
            Zvars_alpha[i] = z_var_alpha;
            z_expr += (1 - z_var_alpha);
            // z_expr += z_var_alpha;

            GRBLinExpr sum_expr_alpha;
            for (int h = 1; h <= ploidy; ++h) {
                if (Zvars_h_alpha[h].find(i) != Zvars_h_alpha[h].end()) {
                    sum_expr_alpha += Zvars_h_alpha[h][i];
                }
            }

            model.addConstr(sum_expr_alpha == ploidy * Zvars_alpha[i], "alpha_constraint_" + std::to_string(i)); // alpha_1 + alpha_2 == 2 * alpha (Pick only one of these)
        }

        for (int32_t i : z_indices_beta) {
            // For beta
            std::string z_var_name_beta = "beta_" + std::to_string(i);
            GRBVar z_var_beta = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, z_var_name_beta);
            Zvars_beta[i] = z_var_beta;
            z_expr += (1 - z_var_beta);
            // z_expr += z_var_beta;

            GRBLinExpr sum_expr_beta;
            for (int h = 1; h <= ploidy; ++h) {
                if (Zvars_h_beta[h].find(i) != Zvars_h_beta[h].end()) {
                    sum_expr_beta += Zvars_h_beta[h][i];
                }
            }

            model.addConstr(sum_expr_beta == Zvars_beta[i], "beta_constraint_" + std::to_string(i)); // beta_1 + beta_2 == beta (Pick only one of these)
        }

        // // Include the z_vars into the objective function
        // for (const auto& z_var_pair : Zvars_alpha) {
        //     z_expr += z_var_pair.second;
        // }
        // for (const auto& z_var_pair : Zvars_beta) {
        //     z_expr += z_var_pair.second;
        // }

        // Print statistics
        fprintf(stderr, "[M::%s::%.3f*%.2f] %.2f%% Minimizers are in ILP\n", __func__,
            realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0),
            (count_kmer_matches * 100.0) / count_sp_r);

        // Clear memory
        Anchor_hits.clear();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Minimizer constraints added to the model\n", __func__,
            realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
        if (is_mixed) {
            fprintf(stderr, "[M::%s::%.3f*%.2f] Using Mixed Integer Programming\n", __func__,
                realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
        } else {
            fprintf(stderr, "[M::%s::%.3f*%.2f] Using Integer Programming\n", __func__,
                realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
        }

        // Continue with the main model construction
        for (int h = 1; h <= ploidy; ++h) {
            // GRBLinExpr vtx_expr;
            // Initialize adjacency lists and expressions for each h
            std::map<std::string, std::vector<std::string>> new_adj;
            GRBLinExpr start_expr;
            GRBLinExpr end_expr;

            // Add start and end variables
            for (int32_t i = 0; i < num_walks; i++) {
                if (paths[i].empty()) continue;
                int32_t u_start = paths[i][0];
                std::string var_name_start = "s_" + std::to_string(u_start) + "_" + std::to_string(i);
                std::string var_name_start_h = var_name_start + "_" + std::to_string(h);

                if (vars[h].find(var_name_start) == vars[h].end()) {
                    GRBVar var_start = model.addVar(0.0, 1.0, 0.0,
                        is_mixed ? GRB_CONTINUOUS : GRB_BINARY, var_name_start_h);
                    vars[h][var_name_start] = var_start;
                }
                start_expr += vars[h][var_name_start];

                int32_t u_end = paths[i].back();
                std::string var_name_end = std::to_string(u_end) + "_" + std::to_string(i) + "_e";
                std::string var_name_end_h = var_name_end + "_" + std::to_string(h);

                if (vars[h].find(var_name_end) == vars[h].end()) {
                    GRBVar var_end = model.addVar(0.0, 1.0, 0.0,
                        is_mixed ? GRB_CONTINUOUS : GRB_BINARY, var_name_end_h);
                    vars[h][var_name_end] = var_end;
                }
                end_expr += vars[h][var_name_end];
            }

            // Set start and end constraints
            model.addConstr(start_expr == 1, "Start_expr_" + std::to_string(h));
            model.addConstr(end_expr == 1, "End_expr_" + std::to_string(h));

            // Without recombination
            for (int32_t i = 0; i < num_walks; i++) {
                if (paths[i].empty()) continue;
                for (int32_t idx = 0; idx < paths[i].size() - 1; idx++) {
                    int32_t u = paths[i][idx];
                    int32_t v = paths[i][idx + 1];
                    std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_"
                        + std::to_string(v) + "_" + std::to_string(i);
                    std::string var_name_h = var_name + "_" + std::to_string(h);

                    // Update adjacency list
                    std::string u_node = std::to_string(u) + "_" + std::to_string(i);
                    std::string v_node = std::to_string(v) + "_" + std::to_string(i);
                    new_adj[u_node].push_back(v_node);

                    if (vars[h].find(var_name) == vars[h].end()) { // Variable does not exist
                        GRBVar var = model.addVar(0.0, 1.0, 0.0,
                            is_mixed ? GRB_CONTINUOUS : GRB_BINARY, var_name_h);
                        vars[h][var_name] = var;
                    }
                    // No need to add 0 * var to the objective
                }
            }

            // Populate hash table to quickly get vertex index in haplotypes
            std::vector<std::unordered_map<int, int>> elementIndexMaps(paths.size());
            for (int h_idx = 0; h_idx < paths.size(); h_idx++) {
                std::unordered_map<int, int> elementIndexMap;
                for (int i = 0; i < paths[h_idx].size(); ++i) {
                    elementIndexMap[paths[h_idx][i]] = i;
                }
                elementIndexMaps[h_idx] = elementIndexMap;
            }

            // Add recombination vertices and edges
            for (int32_t u = 0; u < adj_list.size(); u++) {
                for (auto v : adj_list[u]) {
                    std::string new_vtx = "w_" + std::to_string(u) + "_" + std::to_string(v);

                    bool new_vertex_used = false;
                    for (auto h_idx : haps[u]) {
                        // Check if the next entry paths[h_idx] after u is v
                        int index = elementIndexMaps[h_idx][u];

                        if (index == paths[h_idx].size() - 1 || paths[h_idx][index + 1] != v) {
                            new_vertex_used = true;

                            std::string var_name_1 = std::to_string(u) + "_" + std::to_string(h_idx) + "_"
                                + new_vtx;
                            std::string var_name_1_h = var_name_1 + "_" + std::to_string(h);

                            std::string u_node = std::to_string(u) + "_" + std::to_string(h_idx);
                            new_adj[u_node].push_back(new_vtx);

                            if (vars[h].find(var_name_1) == vars[h].end()) {
                                GRBVar var = model.addVar(0.0, 1.0, 0.0,
                                    is_mixed ? GRB_CONTINUOUS : GRB_BINARY, var_name_1_h);
                                vars[h][var_name_1] = var;
                            }
                            vtx_expr += (c_1 / 2) * vars[h][var_name_1];
                            // vtx_expr += vars[h][var_name_1];
                        }
                    }

                    if (new_vertex_used) {
                        for (auto h_idx : haps[v]) {
                            std::string var_name_2 = new_vtx + "_" + std::to_string(v) + "_"
                                + std::to_string(h_idx);
                            std::string var_name_2_h = var_name_2 + "_" + std::to_string(h);

                            new_adj[new_vtx].push_back(std::to_string(v) + "_" + std::to_string(h_idx));

                            if (vars[h].find(var_name_2) == vars[h].end()) {
                                GRBVar var = model.addVar(0.0, 1.0, 0.0,
                                    is_mixed ? GRB_CONTINUOUS : GRB_BINARY, var_name_2_h);
                                vars[h][var_name_2] = var;
                            }
                            vtx_expr += (c_1 / 2) * vars[h][var_name_2];
                            // vtx_expr += vars[h][var_name_2];
                        }
                    }
                }
            }
            elementIndexMaps.clear();

            // Create the reverse adjacency list
            std::map<std::string, std::vector<std::string>> in_nodes_new;
            for (auto& v : new_adj) {
                for (auto& u : v.second) {
                    in_nodes_new[u].push_back(v.first);
                }
            }

            // Paths-based flow constraints
            for (int32_t i = 0; i < num_walks; i++) {
                if (paths[i].empty()) continue;
                for (int32_t idx = 0; idx < paths[i].size(); idx++) {
                    if (idx == 0 || idx == paths[i].size() - 1) continue; // Skip source and sink nodes

                    GRBLinExpr in_expr;
                    GRBLinExpr out_expr;

                    int32_t v = paths[i][idx];
                    std::string vtx = std::to_string(v) + "_" + std::to_string(i);

                    for (auto& u : in_nodes_new[vtx]) {
                        std::string edge = u + "_" + vtx;
                        if (vars[h].find(edge) != vars[h].end())
                            in_expr += vars[h][edge];
                    }
                    for (auto& u : new_adj[vtx]) {
                        std::string edge = vtx + "_" + u;
                        if (vars[h].find(edge) != vars[h].end())
                            out_expr += vars[h][edge];
                    }

                    std::string constraint_name = "Flow_conservation_" + std::to_string(v) + "_"
                        + std::to_string(i) + "_" + std::to_string(h);
                    model.addConstr(in_expr == out_expr, constraint_name);
                }
            }

            // Recombination vertex flow constraints
            for (int32_t u = 0; u < n_vtx; u++) {
                for (auto v : adj_list[u]) {
                    // For w_u_v vertices
                    std::string w_vtx = "w_" + std::to_string(u) + "_" + std::to_string(v);
                    if (new_adj.find(w_vtx) != new_adj.end()) { // w_vtx exists
                        GRBLinExpr in_expr;
                        GRBLinExpr out_expr;

                        for (auto& adj_node : new_adj[w_vtx]) {
                            std::string edge = w_vtx + "_" + adj_node;
                            if (vars[h].find(edge) != vars[h].end())
                                out_expr += vars[h][edge];
                        }
                        for (auto& adj_node : in_nodes_new[w_vtx]) {
                            std::string edge = adj_node + "_" + w_vtx;
                            if (vars[h].find(edge) != vars[h].end())
                                in_expr += vars[h][edge];
                        }
                        std::string constraint_name = "Flow_conservation_" + w_vtx + "_" + std::to_string(h);
                        model.addConstr(in_expr == out_expr, constraint_name);
                    }
                }
            }

            // Flow constraints for source nodes
            for (int32_t i = 0; i < num_walks; i++) {
                if (paths[i].empty()) continue;
                int32_t u = paths[i][0];
                GRBLinExpr s_expr;
                std::string var_name_start = "s_" + std::to_string(u) + "_" + std::to_string(i);
                s_expr += vars[h][var_name_start];

                std::string vtx = std::to_string(u) + "_" + std::to_string(i);
                for (auto& v : new_adj[vtx]) {
                    std::string edge = vtx + "_" + v;
                    if (vars[h].find(edge) != vars[h].end())
                        s_expr -= vars[h][edge];
                }
                std::string constraint_name = "Source_conservation_" + std::to_string(u) + "_"
                    + std::to_string(i) + "_" + std::to_string(h);
                model.addConstr(s_expr == 0, constraint_name);
            }

            // Flow constraints for sink nodes
            for (int32_t i = 0; i < num_walks; i++) {
                if (paths[i].empty()) continue;
                int32_t u = paths[i].back();
                GRBLinExpr e_expr;
                std::string var_name_end = std::to_string(u) + "_" + std::to_string(i) + "_e";
                std::string vtx = std::to_string(u) + "_" + std::to_string(i);
                for (auto& v : in_nodes_new[vtx]) {
                    std::string edge = v + "_" + vtx;
                    if (vars[h].find(edge) != vars[h].end())
                        e_expr += vars[h][edge];
                }
                e_expr -= vars[h][var_name_end];
                std::string constraint_name = "Sink_conservation_" + std::to_string(u) + "_"
                    + std::to_string(i) + "_" + std::to_string(h);
                model.addConstr(e_expr == 0, constraint_name);
            }

            // Recombination vertex expression <= R;
            // int32_t R = recombination; // Set the maximum number of recombination vertices allowed
            // model.addConstr(vtx_expr <= 2 * recombination, "Max_recombination_constraint" + std::to_string(h)); // there are two edges for 1 recombination
            // Clean up
            new_adj.clear();
            in_nodes_new.clear();
        }
        

        // Set the objective function
        obj = vtx_expr + z_expr;
        model.setObjective(obj, GRB_MINIMIZE);

        // // Maximization objective
        // obj = z_expr;
        // model.setObjective(obj, GRB_MAXIMIZE);

        // Clear data structures if needed
        vars.clear();
        Zvars_h_alpha.clear();
        Zvars_alpha.clear();
        Zvars_h_beta.clear();
        Zvars_beta.clear();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Optimized expanded graph constructed\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));


        // Check the default optimality tolerance
        double defaultTol = model.get(GRB_DoubleParam_OptimalityTol);
        std::cout << "Default Optimality Tolerance: " << defaultTol << std::endl;
        // Set optimality tolerance
        model.set(GRB_DoubleParam_OptimalityTol, 1.00e-08);
        // Optimize model
        model.optimize();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Model optimized\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        // Print constraints
        if (debug)
        {
            printObjectiveFunction(model);
            printConstraints(model);
            printQuadraticConstraints(model);
            printNonZeroVariables(model);
        }

        // Loop over haplotypes h = {1, 2}
        for (int h = 1; h <= ploidy; ++h) {
            // Vector to store the names of non-zero binary variables for haplotype h
            std::vector<std::string> path_strs;

            // Get the list of variables in the model
            GRBVar* variables = model.getVars();
            int num_vars = model.get(GRB_IntAttr_NumVars);

            for (int i = 0; i < num_vars; ++i) {
                GRBVar var = variables[i];
                // Check if the variable is binary and non-zero
                if (var.get(GRB_DoubleAttr_X) == 1.0 || var.get(GRB_DoubleAttr_X) == 1) {
                    std::string var_name = var.get(GRB_StringAttr_VarName);

                    // Skip variables that start with 's', 'a', 'b', or end with 'e'
                    if (var_name[0] == 's' || var_name[0] == 'a' || var_name[0] == 'b' || var_name[var_name.size() - 3] == 'e') {
                        continue;
                    }
                    
                    // std::cerr << "Variable: " << var_name << " = " << var.get(GRB_DoubleAttr_X) << std::endl;
                    if (std::stoi(std::string(1, var_name[var_name.size() - 1])) == h)
                    {
                        std::string var_name_wo_hap = var_name.substr(0, var_name.size() - 2);
                        path_strs.push_back(var_name_wo_hap);
                        // std::cerr << "var_name : " << var_name << "Path string: " << var_name_wo_hap << " On hap: " << std::stoi(std::string(1, var_name[var_name.size() - 1])) << std::endl;
                    }
                }
            }

            // Print path strings
            std::set<std::pair<int32_t, int32_t>> path_vertices_hap;
            int32_t recombination_count = 0;
            for (int i = 0; i < path_strs.size(); i++)
            {
                std::stringstream ss(path_strs[i]);
                std::vector<std::string> tokens;
                std::string item;
                while (std::getline(ss, item, '_')) {
                    tokens.push_back(item);
                }

                std::string hap_1;
                std::string hap_2;
                int32_t u;
                int32_t v;

                if (tokens.size() == 4)
                {
                    u = std::stoi(tokens[0]);
                    hap_1 = tokens[1];
                    v = std::stoi(tokens[2]);
                    hap_2 = tokens[3];
                    path_vertices_hap.insert({u, std::stoi(hap_1)});
                    path_vertices_hap.insert({v, std::stoi(hap_2)});
                    if (debug) std::cerr << "(vtx, hap) => " << "(" << u << "," << hap_1 << ")" << std::endl;
                    if (debug) std::cerr << "(vtx, hap) => " << "(" << v << "," << hap_2 << ")" << std::endl;
                } else {
                    if (tokens[2] == "w")
                    {
                        u = std::stoi(tokens[0]);
                        hap_1 = tokens[1];
                        path_vertices_hap.insert({u, std::stoi(hap_1)});
                        if (debug) std::cerr << "(vtx, hap) => " << "(" << u << "," << hap_1 << ")" << std::endl;
                    } else {
                        v = std::stoi(tokens[3]);
                        hap_2 = tokens[4];
                        path_vertices_hap.insert({v, std::stoi(hap_2)});
                        if (debug) std::cerr << "(vtx, hap) => " << "(" << v << "," << hap_2 << ")" << std::endl;
                    }
                }
            }

            std::vector<std::pair<int32_t, int32_t>> path_vertices_hap_vec(path_vertices_hap.begin(), path_vertices_hap.end());
            

            // Sort path_vertices_hap_vec based on topological order
            std::sort(path_vertices_hap_vec.begin(), path_vertices_hap_vec.end(), [&](std::pair<int32_t, int32_t> a, std::pair<int32_t, int32_t> b) {
                return top_order_map[a.first] < top_order_map[b.first];
            });

            int32_t prev_hap = path_vertices_hap_vec[0].second;
            int32_t prev_str_id = 0;
            int32_t str_id = 0;
            std::vector<std::string> hap_st_en_vec;
            str_id += node_seq[path_vertices_hap_vec[0].first].size();
            for (int32_t i = 1; i < path_vertices_hap_vec.size(); ++i)
            {
                str_id += node_seq[path_vertices_hap_vec[i].first].size();

                if (prev_hap != path_vertices_hap_vec[i].second)
                {
                    recombination_count++;
                    std::string str = ">(" + hap_id2name[prev_hap] + ",[" + std::to_string(prev_str_id) + "," + std::to_string(str_id - 1) + "])";
                    hap_st_en_vec.push_back(str);
                    prev_hap = path_vertices_hap_vec[i].second;
                    prev_str_id = str_id;
                }
            }

            // Capture the last segment after the loop
            std::string str = ">(" + hap_id2name[path_vertices_hap_vec.back().second] + ",[" + std::to_string(prev_str_id) + "," + std::to_string(str_id - 1) + "])";
            hap_st_en_vec.push_back(str);

            std::cerr << "Recombination count for haplotype " << h << ": " << recombination_count << std::endl;
            if (recombination_count > 0)
            {
                std::cerr << "Recombined haplotypes for haplotype " << h << ": ";
                for (int i = 0; i < hap_st_en_vec.size(); i++)
                {
                    std::cerr << hap_st_en_vec[i];
                }
                std::cerr << std::endl;
            } else {
                std::cerr << "Recombined haplotypes for haplotype " << h << ": ";
                int32_t sum_str = 0;
                for (int i = 0; i < path_vertices_hap_vec.size(); i++)
                {
                    sum_str += node_seq[path_vertices_hap_vec[i].first].size();
                }
                std::cerr << ">(" << hap_id2name[prev_hap] << ",[" << 0 << "," << sum_str - 1 << "])" << std::endl;
            }

            // Verify the path vertices by checking if there exists an edge between the vertices
            if (debug) std::cout << "(" << "s" << "," << path_vertices_hap_vec[0].first << ")" << "->";
            for (int i = 1; i < path_vertices_hap_vec.size(); i++)
            {
                int32_t u = path_vertices_hap_vec[i - 1].first;
                int32_t v = path_vertices_hap_vec[i].first;
                bool exist_edge = false;
                for (auto w : adj_list[u])
                {
                    if (w == v)
                    {
                        exist_edge = true;
                        break;
                    }
                }
                if (!exist_edge)
                {
                    fprintf(stderr, "Error: No edge between %d and %d\n", u, v);
                    exit(1);
                }
                if (debug) std::cout << "(" << u << "," << v << ")" << "->";
            }
            if (debug) std::cout << "(" << path_vertices_hap_vec.back().first << "," << "e" << ")" << std::endl;

            // Get the path string and store in haplotype[h]
            for (int i = 0; i < path_vertices_hap_vec.size(); i++)
            {
                haplotypes[h] += node_seq[path_vertices_hap_vec[i].first];
            }
        } // End of loop over h

    } catch (GRBException e) {
        std::cerr << "Error code = " << e.getErrorCode() << std::endl;
        std::cerr << e.getMessage() << std::endl;
    } catch (...) {
        std::cerr << "Exception during optimization" << std::endl;
    }

    // After the try-catch block, write the haplotypes to files
    for (int h = 1; h <= ploidy; ++h) {
        // Write haplotype to a file as fasta from the path
        std::string path_str = haplotypes[h];
        std::string hap_file_ = hap_file + "_" + std::to_string(h) + ".fa";
        std::string hap_name_ = hap_name + "_" + std::to_string(h);
        std::ofstream hap_file_stream(hap_file_, std::ios::out);
        hap_file_stream << ">" << hap_name_ << " LN:" << path_str.size() << std::endl;
        // Write the path_str to the file 80 characters per line
        for (size_t i = 0; i < path_str.size(); i += 80) {
            hap_file_stream << path_str.substr(i, 80) << std::endl;
        }
        hap_file_stream.close();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Haplotype %d of size: %lu written to: %s\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), h, path_str.size(), hap_file_.c_str());
    }
}
