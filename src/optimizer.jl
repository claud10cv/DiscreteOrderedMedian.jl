function glpk_optimizer_data(opt)
    attrs = (m, t, _) -> set_optimizer_attributes(m, "msg_lev" => GLPK.GLP_MSG_OFF, "tm_lim" => t)
    return OptimizerData(opt, attrs)
end

function cplex_optimizer_data(opt)
    attrs = (m, t, sl) -> set_optimizer_attributes(m, 
                                            "CPXPARAM_ScreenOutput" => 0,
                                            "CPXPARAM_TimeLimit" => t,
                                            "CPXPARAM_Threads" => 1,
                                            "CPXPARAM_MIP_Tolerances_MIPGap" => 0.0,
                                            "CPXPARAM_MIP_Limits_Solutions" => sl)
    return OptimizerData(opt, attrs)
end

function gurobi_optimizer_data(opt)
    attrs = (m, t, sl) -> set_optimizer_attributes(m, 
                                            "OutputFlag" => 0,
                                            "Threads" => 1,
                                            "TimeLimit" => t,
                                            "MIPGap" => 0.0,
                                            "SolutionLimit" => sl)
    return OptimizerData(opt, attrs)
end

default_parameters() = Parameters(true, true, true, true, 7200.0, glpk_optimizer_data(GLPK.Optimizer))