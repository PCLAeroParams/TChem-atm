for number_of_particles in ${nparticles[@]}; do
for N in ${nbatch[@]}; do
  for threads in ${Nthread[@]}; do
      echo "nbatch = $N"
      echo "Nthread = $threads"
      thread_size=${threads}
      output_wall_times="${tchem_outputs}/${experiment_name}/wall_times_nbatch_${N}_thread_size_${thread_size}_number_of_particles${number_of_particles}.json"
      output="${tchem_outputs}/${experiment_name}/reaction_rates_nbatch_${N}_thread_size_${thread_size}.txt"
      echo "${tchem_outputs}"
      echo "${experiment_name}"
      echo "/reaction_rates_nbatch_${N}_thread_size_${thread_size}.txt"
      echo "output file = " $output
      run_this="OMP_NUM_THREADS=$thread_size OMP_PLACES=threads OMP_PROC_BIND=close ${exec} \
               --batch_size=$N \
               --use_cloned_samples=$use_cloned_samples \
               --verbose=$verbose \
               --outputfile_times=$output_wall_times \
               --outputfile=$output \
               --number_of_particles=$number_of_particles \
               $scenario_n_inputs \
               $numerical_params "
      echo $run_this
      eval $run_this
      sleep 2
      # Note: we use this with kokkos tools
      kp_json_writer $machine_name* > "${tchem_outputs}/${experiment_name}/simple_timer_nbatch_${N}_thread_size_${thread_size}_number_of_particles${number_of_particles}.json"
      sleep 2
      rm -rf $machine_name*
    done
done
done