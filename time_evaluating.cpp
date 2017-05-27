////
//// Created by Olia on 01.05.2017.
////
//
//#include <chrono>
//#include <atomic>
//using namespace std;
//
//inline std::chrono::high_resolution_clock::time_point get_current_time_fenced()
//{
//    std::atomic_thread_fence(memory_order_seq_cst);
//    auto res_time = std::chrono::high_resolution_clock::now();
//    std::atomic_thread_fence(memory_order_seq_cst);
//    return res_time;
//}
//
//template<class D>
//inline long long to_us(const D& d)
//{
//    return std::chrono::duration_cast<chrono::microseconds>(d).count();
//}