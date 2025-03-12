for N in ${nbatch[@]}; do
  for i in "${!vector_size[@]}"; do
    for j in "${!team_size[@]}"; do
      vector_thread_size=${vector_size[i]}
      team_thread_size=${team_size[j]}
      #if ((vector_thread_size*$team_thread_size < 1024)); then
        echo "nbatch = $N"
        echo "vector thread size = $vector_thread_size"
        echo "team thread size = $team_thread_size"
        output_wall_times="${tchem_outputs}/${experiment_name}/wall_times_nbatch_${N}_vecsize_${vector_thread_size}_teamThread_size_${team_thread_size}.json"
        output_file="${tchem_outputs}/${experiment_name}/reaction_rates_nbatch_${N}_vecsize_${vector_thread_size}_teamThread_size_${team_thread_size}.txt"
        echo "${tchem_outputs}"
        echo "${experiment_name}"
        echo "/reaction_rates_nbatch_${N}_vecsize_${vector_thread_size}_teamThread_size_${team_thread_size}.txt"
        echo "output file = ${output_file}"
        run_this="OMP_NUM_THREADS=16 OMP_PLACES=threads OMP_PROC_BIND=close $exec \
                --batch_size=$N \
                --team_thread_size=$team_thread_size \
                --vector_thread_size=$vector_thread_size \
                --use_cloned_samples=$use_cloned_samples \
                --verbose=$verbose \
                --outputfile_times=$output_wall_times \
                --outputfile=$output_file \
                $scenario_n_inputs \
                $numerical_params "
        echo $run_this
        eval $run_this
        sleep 2
        kp_json_writer $machine_name* > "${tchem_outputs}/${experiment_name}/simple_timer_nbatch_${N}_vecsize_${vector_thread_size}_teamThread_size_${team_thread_size}.json"
        sleep 2
        rm -rf $machine_name*
      #fi
done
done
done

