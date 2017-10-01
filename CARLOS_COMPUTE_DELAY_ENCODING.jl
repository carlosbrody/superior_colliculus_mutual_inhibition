end_of_delay = Int(round(F["rule_and_delay_periods"]/model_params[:dt]))
pronodes =  proVall[[1,4],end_of_delay,:];
antinodes = proVall[[2,3],end_of_delay,:];

mean_pro_encode = round(mean(pronodes)*1000)/1000;
mean_anti_encode = round(mean(antinodes)*1000)/1000;
std_pro_encode = round(std(pronodes)*1000)/1000;
std_anti_encode = round(std(antinodes)*1000)/1000;


print("Average PRO  node external activity at end of delay period was ");print(mean_pro_encode);print(" w/ std +/- ");print(std_pro_encode);   println()
print("Average ANTI node external activity at end of delay period was ");print(mean_anti_encode);print(" w/ std +/- ");print(std_anti_encode);   println()

