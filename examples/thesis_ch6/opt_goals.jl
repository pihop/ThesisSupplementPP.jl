using Distributions
# What does the reward function aim to accomplish? 
# 80% of the swarm at the target location. Minimise the sensing rate? 

function goal_is_sat(traj, var)
    return maximum(map(x -> (x >= 0.8), traj))
end

function prob_goal_is_sat(traj, var)
    dist = map((m, s) -> Normal(m, sqrt(s)), traj, var) 
    return maximum(map(x -> (1-cdf(x, 0.8)) > 0.9, dist))
end
