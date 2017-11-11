/*  Sunil
 *  Umasankar
 *  suniluma
 */

#ifndef A1_HPP
#define A1_HPP

#include <vector>
#include <set>
#include <mpi.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iterator>
using namespace std;

int connected_components(std::vector<signed char>& A, int n, int q, const char* out, MPI_Comm comm) {
    // ...
    int rank;
    MPI_Comm_rank(comm, &rank);

    int max = n/q;
    int vec_len = max*max;
    vector<int> P(vec_len,0);
    vector<int> Pl(vec_len,0);
    vector<int> M(vec_len,0);
    vector<int> Q(vec_len,0);

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if(A[i*max+j] == 1) {
                P[i*max+j] = (rank/q)*max + i;
            }
        }
    }


    vector <int> high(max,0);


    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if (P[i + j*max] > high[i]) {
                high[i] = P[i + j*max];
            }
        }
    }


    int col = rank % q;
    vector<int> back_high(max,0);
    MPI_Comm col_world;
    MPI_Comm_split(comm, col, rank, &col_world);
    MPI_Allreduce(high.data(),back_high.data(),max,MPI_INT,MPI_MAX,col_world);

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
                P[i*max+j] = back_high[j];
        }
    }

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if(A[i*max+j] == 1) {
                M[i*max+j] = P[i*max+j];
            }
        }
    }

/*  Testing:

    if(rank == 8) {
        for(int i = 0; i < max; i++) {
            for (int j = 0; j < max; j++) {
                cout<<M[i*max+j];
            }
            cout<<"\n";
        }
    }
*/

    vector <int> high_row(max,0);

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if (M[j + i*max] > high_row[i]) {
                high_row[i] = M[j + i*max];
            }
        }
    }

    int row = rank / q;
    vector <int> back_high_row(max,0);
    MPI_Comm row_world;
    MPI_Comm_split(comm, row, rank, &row_world);
    MPI_Allreduce(high_row.data(),back_high_row.data(),max,MPI_INT,MPI_MAX,row_world);




    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            Q[i*max+j] = back_high_row[i];
            if(Q[i*max+j] == (rank%q)*max+j) {
                M[i*max+j] = P[i*max+j];
            }
        }
    }

    vector <int> p_prime_local(max,0);

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if (M[j + i*max] > p_prime_local[i]) {
                p_prime_local[i] = M[j + i*max];
            }
        }
    }

    vector <int> p_prime_global(max,0);
    MPI_Comm prime_world;
    MPI_Comm_split(comm, row, rank, &prime_world);
    MPI_Allreduce(p_prime_local.data(),p_prime_global.data(),max,MPI_INT,MPI_MAX,prime_world);


    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            Pl[i*max+j] = p_prime_global[i];
        }
    }

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if(P[i*max+j] == (rank/q)*max + i) {
                M[i*max+j] = Pl[i*max+j];
            }
        }
    }


    vector <int> m_local(max,0);

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if (M[j + i*max] > m_local[i]) {
                m_local[i] = M[j + i*max];
            }
        }
    }

    vector <int> m_global(max,0);
    MPI_Comm m_world;
    MPI_Comm_split(comm, row, rank, &m_world);
    MPI_Allreduce(m_local.data(),m_global.data(),max,MPI_INT,MPI_MAX,m_world);

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            Q[i*max+j] = m_global[i];
        }
    }

    for(int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if(Pl[i*max+j] > Q[i*max+j]) {
                P[i*max+j] = Pl[i*max+j];
            } else {
                P[i*max+j] = Q[i*max+j];
            }
        }
    }

    vector <int> all_size(n*n,0);
    MPI_Allgather(P.data(),max*max,MPI_INT,all_size.data(),max*max,MPI_INT,comm);
    set<int> s(all_size.begin(),all_size.end());
    if(rank == 0) cout<<s.size();

    vector <int> output;
    for(int i = 0; i < n*n; i++) {
        if (i%n <= q) {
            if(i%q == 0) {
                output.push_back(all_size[i]);
            }
        }
    }

    if(rank == 0) {
        std::ofstream FILE(out, std::ios::out | std::ofstream::binary);
        std::copy(output.begin(), output.end(), std::ostreambuf_iterator<char>(FILE));
    }

    return s.size();
} // connected_components

#endif // A1_HPP
