function step_through(d::ncData)
    for t in 1:length(d.timevec):
        mean(d.data[:,:,t])
    end
end