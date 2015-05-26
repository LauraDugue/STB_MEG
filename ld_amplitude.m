function ld_amplitude(dirmat,obs,cond)

for elec = 1:64

    dirElec = [dirmat '/timeFreq/' obs '_elec' num2str(elec) '_' cond{:}];
    load(dirElec)
    amplitude = mean(abs(tf),3);
    save([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_' cond{:} '.mat'],'amplitude','freqs','times');
    
end

end