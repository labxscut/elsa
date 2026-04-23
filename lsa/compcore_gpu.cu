#include "compcore.hpp"
#include <cuda_runtime.h>
#include <climits>
#include <stdexcept>

using namespace std;

namespace {

__global__ void lla_row_kernel(const double* x,
                               const double* y,
                               const double* z,
                               int n,
                               int max_shift,
                               int row_i,
                               const double* prev_psm,
                               const double* prev_nsm,
                               const int* prev_lpsm,
                               const int* prev_lnsm,
                               double* curr_psm,
                               double* curr_nsm,
                               int* curr_lpsm,
                               int* curr_lnsm) {
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if (k > n) {
        return;
    }

    if (max_shift != INT_MAX && abs(row_i - k) > max_shift) {
        curr_psm[k] = 0.0;
        curr_nsm[k] = 0.0;
        curr_lpsm[k] = 0;
        curr_lnsm[k] = 0;
        return;
    }

    double s1 = x[row_i - 1] * y[row_i - 1] * z[k - 1];

    double prev_p = (k > 1) ? prev_psm[k - 1] : 0.0;
    double prev_n = (k > 1) ? prev_nsm[k - 1] : 0.0;
    int prev_lp = (k > 1) ? prev_lpsm[k - 1] : 0;
    int prev_ln = (k > 1) ? prev_lnsm[k - 1] : 0;

    double cand_p = prev_p + s1;
    if (cand_p > 0.0) {
        curr_psm[k] = cand_p;
        curr_lpsm[k] = prev_lp + 1;
    } else {
        curr_psm[k] = 0.0;
        curr_lpsm[k] = 0;
    }

    double cand_n = prev_n - s1;
    if (cand_n > 0.0) {
        curr_nsm[k] = cand_n;
        curr_lnsm[k] = prev_ln + 1;
    } else {
        curr_nsm[k] = 0.0;
        curr_lnsm[k] = 0;
    }
}

template <typename T>
void cuda_copy_to_device(T*& device_ptr, const T* host_ptr, size_t count) {
    cudaMalloc(reinterpret_cast<void**>(&device_ptr), count * sizeof(T));
    cudaMemcpy(device_ptr, host_ptr, count * sizeof(T), cudaMemcpyHostToDevice);
}

template <typename T>
void cuda_allocate_zeroed(T*& device_ptr, size_t count) {
    cudaMalloc(reinterpret_cast<void**>(&device_ptr), count * sizeof(T));
    cudaMemset(device_ptr, 0, count * sizeof(T));
}

bool is_better_candidate(double candidate_abs,
                         int candidate_len,
                         int candidate_i,
                         int candidate_j,
                         int candidate_k,
                         double best_abs,
                         int best_len,
                         int best_i,
                         int best_j,
                         int best_k) {
    if (candidate_abs > best_abs) {
        return true;
    }

    if (candidate_abs == best_abs) {
        if (candidate_len > 0 && (best_len == 0 || candidate_len > best_len)) {
            return true;
        }

        if (candidate_len == best_len && best_len > 0) {
            if (candidate_i < best_i) {
                return true;
            }
            if (candidate_i == best_i && candidate_j < best_j) {
                return true;
            }
            if (candidate_i == best_i && candidate_j == best_j && candidate_k < best_k) {
                return true;
            }
        }
    }

    return false;
}

}

int test() {
    cout << "tested";
    return 0;
}

double calc_LA(VectorDouble x, VectorDouble y, VectorDouble z) {
    if (x.size() != y.size() || y.size() != z.size()) {
        throw std::runtime_error("All input vectors must have the same length");
    }

    double sum = 0.0;
    int n = static_cast<int>(x.size());
    for (int i = 0; i < n; i++) {
        sum += (x[i] * y[i] * z[i]);
    }
    return sum / n;
}

LSA_Result DP_lsa(const LSA_Data& data, bool keep_trace) {
    LSA_Result lsa_result;
    int max_p[2] = {0};
    int porn = 0;
    double max_s = -std::numeric_limits<double>::infinity();
    MatrixDouble psm = vector<vector<double>>(data.X.size() + 1, vector<double>(data.Y.size() + 1));
    MatrixDouble nsm = vector<vector<double>>(data.X.size() + 1, vector<double>(data.Y.size() + 1));

    for (unsigned int i = 1; i <= data.X.size(); i++) {
        for (unsigned int j = std::max(1, (int)i - data.max_shift); (int)j <= std::min((int)data.Y.size(), (int)i + data.max_shift); j++) {
            double s1 = data.X[i - 1] * data.Y[j - 1];
            psm[i][j] = std::max(0., psm[i - 1][j - 1] + s1);
            nsm[i][j] = std::max(0., nsm[i - 1][j - 1] - s1);
            if (psm[i][j] >= max_s) {
                max_p[0] = i;
                max_p[1] = j;
                max_s = psm[i][j];
                porn = 1;
            }
            if (nsm[i][j] >= max_s) {
                max_p[0] = i;
                max_p[1] = j;
                max_s = nsm[i][j];
                porn = -1;
            }
        }
    }

    int length = 0;
    vector<int> step;
    step.push_back(max_p[0]);
    step.push_back(max_p[1]);
    if (porn == -1) {
        lsa_result.score = -1 * nsm[max_p[0]][max_p[1]] / data.X.size();
        while (nsm[max_p[0] - length][max_p[1] - length] != 0. && keep_trace == true) {
            length++;
            lsa_result.trace.push_back(step);
            step.clear();
            step.push_back(max_p[0] - length);
            step.push_back(max_p[1] - length);
        }
    } else {
        lsa_result.score = psm[max_p[0]][max_p[1]] / data.X.size();
        while (psm[max_p[0] - length][max_p[1] - length] != 0. && keep_trace == true) {
            length++;
            lsa_result.trace.push_back(step);
            step.clear();
            step.push_back(max_p[0] - length);
            step.push_back(max_p[1] - length);
        }
    }

    return lsa_result;
}

LLA_Result DP_lla(const LLA_Data& data, bool keep_trace) {
    LLA_Result lla_result;

    if (data.X.empty() || data.Y.empty() || data.Z.empty()) {
        lla_result.score = 0.0;
        return lla_result;
    }
    if (data.X.size() != data.Y.size() || data.Y.size() != data.Z.size()) {
        lla_result.score = 0.0;
        return lla_result;
    }

    const int n = static_cast<int>(data.X.size());
    const size_t vec_count = static_cast<size_t>(n + 1);

    double* d_x = nullptr;
    double* d_y = nullptr;
    double* d_z = nullptr;
    double* d_prev_psm = nullptr;
    double* d_prev_nsm = nullptr;
    double* d_curr_psm = nullptr;
    double* d_curr_nsm = nullptr;
    int* d_prev_lpsm = nullptr;
    int* d_prev_lnsm = nullptr;
    int* d_curr_lpsm = nullptr;
    int* d_curr_lnsm = nullptr;

    cuda_copy_to_device(d_x, data.X.data(), data.X.size());
    cuda_copy_to_device(d_y, data.Y.data(), data.Y.size());
    cuda_copy_to_device(d_z, data.Z.data(), data.Z.size());
    cuda_allocate_zeroed(d_prev_psm, vec_count);
    cuda_allocate_zeroed(d_prev_nsm, vec_count);
    cuda_allocate_zeroed(d_curr_psm, vec_count);
    cuda_allocate_zeroed(d_curr_nsm, vec_count);
    cuda_allocate_zeroed(d_prev_lpsm, vec_count);
    cuda_allocate_zeroed(d_prev_lnsm, vec_count);
    cuda_allocate_zeroed(d_curr_lpsm, vec_count);
    cuda_allocate_zeroed(d_curr_lnsm, vec_count);

    vector<double> h_curr_psm(vec_count, 0.0);
    vector<double> h_curr_nsm(vec_count, 0.0);
    vector<int> h_curr_lpsm(vec_count, 0);
    vector<int> h_curr_lnsm(vec_count, 0);

    double best_abs = 0.0;
    int best_len = 0;
    int best_sign = 1;
    int best_i = 0;
    int best_j = 0;
    int best_k = 0;

    constexpr int block_size = 256;
    dim3 block(block_size);
    dim3 grid((n + block_size - 1) / block_size);

    for (int row_i = 1; row_i <= n; ++row_i) {
        lla_row_kernel<<<grid, block>>>(d_x,
                                        d_y,
                                        d_z,
                                        n,
                                        data.max_shift,
                                        row_i,
                                        d_prev_psm,
                                        d_prev_nsm,
                                        d_prev_lpsm,
                                        d_prev_lnsm,
                                        d_curr_psm,
                                        d_curr_nsm,
                                        d_curr_lpsm,
                                        d_curr_lnsm);
        cudaDeviceSynchronize();

        cudaMemcpy(h_curr_psm.data(), d_curr_psm, vec_count * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_curr_nsm.data(), d_curr_nsm, vec_count * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_curr_lpsm.data(), d_curr_lpsm, vec_count * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_curr_lnsm.data(), d_curr_lnsm, vec_count * sizeof(int), cudaMemcpyDeviceToHost);

        for (int k = 1; k <= n; ++k) {
            double current_abs = 0.0;
            int current_len = 0;
            int current_sign = 1;

            if (h_curr_psm[k] >= h_curr_nsm[k]) {
                current_abs = h_curr_psm[k];
                current_len = h_curr_lpsm[k];
                current_sign = 1;
            } else {
                current_abs = h_curr_nsm[k];
                current_len = h_curr_lnsm[k];
                current_sign = -1;
            }

            if (current_abs > 0.0 && is_better_candidate(current_abs, current_len, row_i, row_i, k,
                                                         best_abs, best_len, best_i, best_j, best_k)) {
                best_abs = current_abs;
                best_len = current_len;
                best_sign = current_sign;
                best_i = row_i;
                best_j = row_i;
                best_k = k;
            }
        }

        swap(d_prev_psm, d_curr_psm);
        swap(d_prev_nsm, d_curr_nsm);
        swap(d_prev_lpsm, d_curr_lpsm);
        swap(d_prev_lnsm, d_curr_lnsm);

        cudaMemset(d_curr_psm, 0, vec_count * sizeof(double));
        cudaMemset(d_curr_nsm, 0, vec_count * sizeof(double));
        cudaMemset(d_curr_lpsm, 0, vec_count * sizeof(int));
        cudaMemset(d_curr_lnsm, 0, vec_count * sizeof(int));
    }

    if (best_abs == 0.0) {
        lla_result.score = 0.0;
    } else {
        lla_result.score = (best_sign == 1 ? best_abs : -best_abs) / data.X.size();
    }

    if (keep_trace && best_len > 0) {
        int start_i = best_i - best_len + 1;
        int start_k = best_k - best_len + 1;
        lla_result.trace.push_back({best_i, best_j, best_k});
        lla_result.trace.push_back({start_i, start_i, start_k});
    }

    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_z);
    cudaFree(d_prev_psm);
    cudaFree(d_prev_nsm);
    cudaFree(d_curr_psm);
    cudaFree(d_curr_nsm);
    cudaFree(d_prev_lpsm);
    cudaFree(d_prev_lnsm);
    cudaFree(d_curr_lpsm);
    cudaFree(d_curr_lnsm);

    return lla_result;
}
