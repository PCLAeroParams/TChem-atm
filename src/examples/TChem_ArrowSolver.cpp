#include <Kokkos_Core.hpp> // technically this is encapsulated by TChem.hpp
#include <KokkosBlas.hpp> // For BLAS functions like gemm
#include "TChem.hpp"

// Function to print a 2D host view
template<class HostViewType>
void print_2d_view(HostViewType hv) {
    // It's a good idea to check the dimensions before printing
    // to avoid flooding the console with enormous amounts of data.
    if (hv.extent(0) > 20 || hv.extent(1) > 20) {
        printf("View is too large to print.\n");
        printf("Dimensions: %zu x %zu\n", hv.extent(0), hv.extent(1));
        return;
    }

    printf("------------------------\n");
    for (size_t i = 0; i < hv.extent(0); ++i) {
        for (size_t j = 0; j < hv.extent(1); ++j) {
            printf("%6.2f ", hv(i, j));
        }
        printf("\n");
    }
    printf("------------------------\n");
}

template<class HostViewType>
void print_1d_view(HostViewType hv) {
    printf("------------------------\n");
    for (size_t i = 0; i < hv.extent(0); ++i) {
        printf("%6.2f ", hv(i));
    }
    printf("\n------------------------\n");
}

int main(int argc, char* argv[]){

    Kokkos::initialize(argc, argv);
    {   

        const bool detail = false;
        auto host_exec_space = TChem::host_exec_space();
        auto exec_space = TChem::exec_space();
        TChem::exec_space().print_configuration(std::cout, detail);
        TChem::host_exec_space().print_configuration(std::cout, detail);

        // allocate 5x5 view
        const int nRows = 5;
        const int nColumns = 5;
        Kokkos::View<float**> P_device("P_device", nRows, nColumns);

        // lambda*v1 for v1 the first eigenvector and corres. eigenval
        const std::vector<float> b_data = { 2.12132034f, 1.06066017f, 1.06066017f, 1.06066017f, 1.06066017f };
        const int size = b_data.size();

        Kokkos::View<float*> b_device("b_device", nRows);

        // temporary, unmanaged host view that wraps the source data.
        // This tells Kokkos how to access the host data without managing its memory.
        Kokkos::View<const float*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            b_data_unmanaged(b_data.data(), size);

        // copy data from unmanaged view to b_device view
        Kokkos::deep_copy(exec_space, b_device, b_data_unmanaged);

        // initialize view with arrow configuration
        Kokkos::parallel_for("set_values", nRows, KOKKOS_LAMBDA(const int i) {
            for (int j = 0; j < nColumns; j++){
                if (i == j){
                    P_device(i, j) = 1.0; 
                }
                if (i == 0){
                    P_device(i, j) = 1.0;
                }
                if (j == 0){
                    P_device(i, j) = 1.0;
                }
            }
        });

        // block matrix components
        auto A_device = Kokkos::subview(P_device, std::make_pair(0, 1), std::make_pair(0, 1));
        auto B_device = Kokkos::subview(P_device, std::make_pair(0, 1), std::make_pair(1, nColumns));
        auto C_device = Kokkos::subview(P_device, std::make_pair(1, nRows), std::make_pair(0, 1));
        auto D_device = Kokkos::subview(P_device, std::make_pair(1, nRows), std::make_pair(1, nColumns));
        
        // b vector partitions into q and r vectors 
        auto q_device = Kokkos::subview(b_device, std::make_pair(0, 1));
        auto r_device = Kokkos::subview(b_device, std::make_pair(1, nRows));

        // compute D inverse
        // TODO: May not want to alter the contents of D which points to P
        auto D_inv_device = D_device; // shallow copy (point to memory of D)
        Kokkos::parallel_for("D inverse", nRows-1, KOKKOS_LAMBDA(const int i) {
            D_inv_device(i, i) = 1.0f / D_inv_device(i, i); 
        });

        // --- Compute Schur Complement: S = A - B * D_inv * C ---

        // Step 1: Allocate intermediate views on the device
        Kokkos::View<float**> temp_B_Dinv_device("temp_B_Dinv", 1, nColumns - 1);
        Kokkos::View<float**> result_prod_device("result_prod", 1, 1);

        // Step 2: Compute temp = B * D_inv
        KokkosBlas::gemm(exec_space, "N", "N", 1.0f, B_device, D_inv_device, 0.0f, temp_B_Dinv_device);

        // Step 3: Compute result = temp * C
        KokkosBlas::gemm(exec_space, "N", "N", 1.0f, temp_B_Dinv_device, C_device, 0.0f, result_prod_device);

        // Step 4: Compute S = A - result. This is an element-wise operation.
        // We overwrite A_device with the result, which is the Schur complement.
        Kokkos::parallel_for("element_wise_sub", 1, KOKKOS_LAMBDA(const int& ) {
            // TODO: May not want to alter the contents of A which points to P
            A_device(0, 0) = A_device(0, 0) - result_prod_device(0, 0);
        });
        

        // Compute q - B * D inverse * r


        // create a copy of P on the host with same memory layout as device
        using device_view_type = decltype(P_device);
        Kokkos::View<
            device_view_type::data_type,
            device_view_type::array_layout, // Use the same layout as the device view
            Kokkos::HostSpace               // But place it on the host
        > P_host("P_host", nRows, nColumns);
        Kokkos::deep_copy(exec_space, P_host, P_device);

        auto A_host = Kokkos::subview(P_host, std::make_pair(0, 1), std::make_pair(0, 1));
        auto B_host = Kokkos::subview(P_host, std::make_pair(0, 1), std::make_pair(1, nColumns));
        auto C_host = Kokkos::subview(P_host, std::make_pair(1, nRows), std::make_pair(0, 1));
        auto D_host = Kokkos::subview(P_host, std::make_pair(1, nRows), std::make_pair(1, nColumns));

        printf("P\n");
        print_2d_view(P_host);
        printf("A\n");
        print_2d_view(A_host);
        printf("B\n");
        print_2d_view(B_host);
        printf("C\n");
        print_2d_view(C_host);
        printf("D\n");
        print_2d_view(D_host);

        /*
        auto b_host = Kokkos::create_mirror_view_and_copy(host_exec_space, b_device);
        auto q_host = Kokkos::subview(b_host, std::make_pair(0, 1));
        auto r_host = Kokkos::subview(b_host, std::make_pair(1, nRows));

        printf("q\n");
        print_1d_view(q_host);
        printf("r\n");
        print_1d_view(r_host);
        */

        printf("End of Kokkos execution\n");
    }
    Kokkos::finalize();

    return 0;
}