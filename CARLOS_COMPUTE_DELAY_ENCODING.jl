end_of_delay = Int(round(F["rule_and_delay_periods"]/model_params[:dt]))
p_pronodes =  proVall[[1,4],end_of_delay,:];
p_antinodes = proVall[[2,3],end_of_delay,:];

p_mean_pro_encode = round(mean(p_pronodes)*1000)/1000;
p_mean_anti_encode = round(mean(p_antinodes)*1000)/1000;
p_std_pro_encode = round(std(p_pronodes)*1000)/1000;
p_std_anti_encode = round(std(p_antinodes)*1000)/1000;
p_rule_strength = round((p_mean_pro_encode - p_mean_anti_encode)*1000)/1000;

println("ON PRO TRIALS")
print("Average PRO  node external activity at end of delay period was ");print(p_mean_pro_encode);print(" w/ std +/- ");print(p_std_pro_encode);   println()
print("Average ANTI node external activity at end of delay period was ");print(p_mean_anti_encode);print(" w/ std +/- ");print(p_std_anti_encode);   println()
print("Average PRO  RULE encoding was "); print(p_rule_strength); println()

a_pronodes =  antiVall[[1,4],end_of_delay,:];
a_antinodes = antiVall[[2,3],end_of_delay,:];
a_mean_pro_encode = round(mean(a_pronodes)*1000)/1000;
a_mean_anti_encode = round(mean(a_antinodes)*1000)/1000;
a_std_pro_encode = round(std(a_pronodes)*1000)/1000;
a_std_anti_encode = round(std(a_antinodes)*1000)/1000;
a_rule_strength = round((a_mean_anti_encode - a_mean_pro_encode)*1000)/1000;

println()
println("ON ANTI TRIALS")
print("Average PRO  node external activity at end of delay period was ");print(a_mean_pro_encode);print(" w/ std +/- ");print(a_std_pro_encode);   println()
print("Average ANTI node external activity at end of delay period was ");print(a_mean_anti_encode);print(" w/ std +/- ");print(a_std_anti_encode);   println()
print("Average ANTI RULE encoding was "); print(a_rule_strength); println()
