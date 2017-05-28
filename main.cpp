#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <mutex>
#include <thread>
#include <map>
#include <sstream>

using namespace std;

mutex mx;

inline std::chrono::high_resolution_clock::time_point get_current_time_fenced()
{
    std::atomic_thread_fence(memory_order_seq_cst);
    auto res_time = std::chrono::high_resolution_clock::now();
    std::atomic_thread_fence(memory_order_seq_cst);
    return res_time;
}

template<class D>
inline long long to_us(const D& d)
{
    return std::chrono::duration_cast<chrono::microseconds>(d).count();
}
//Write my func
double func_calculation(double x1, double x2, int m) {
//    double sum1 = 0;
//    double sum2 = 0;
//    double g;
//    for (int i = 1; i <= m; ++i)
//    {
//        sum1 += i * cos((i + 1) * x1 + 1);
//        sum2 += i * cos((i + 1) * x2 + 1);
//    }
//
//    g = - sum1 * sum2;
//
//    return g;
    double g,sum = 0.0;
    int j;
    for (int i = -2; i <= 2; ++i)
    {
        j = i;
        sum += 1 / (5 * (i + 2) + j + 3 + pow(x1 - 16* j,6) + pow(x2 - 16* i,6));
    }

    g = pow(0.002 + sum, -1);

    return g;

}

double integration(double x0, double x, double y0, double y, int m, double pr) {
    assert (m >= 5);
    double sum = 0;
    for (double i = x0; i <= x; i += pr) for (double j = y0; j <= y; j += pr)
            sum += func_calculation(i + pr / 2.0, j + pr / 2.0, m) * pr * pr;
    return sum;
}

void integrateWithThreads(double x0, double x, double y0, double y, int m, double pr, double *r) {
    double integAnswer = integration(x0, x, y0, y, m, pr);
    lock_guard<mutex> lg(mx);
    *r += integAnswer;
}

template <class T>
T get_param(string key, map<string, string> myMap) {
    istringstream ss(myMap[key]);
    T val;
    ss >> val;
    return val;
}

map<string, string> read_config(string filename) {
    string line, delimiter;
    ifstream myfile (filename);
    map<string, string> mp;
    delimiter = "=";
    if (myfile.is_open())
    {
        while (getline(myfile,line))
        {
            int pos = line.find(delimiter);
            string key = line.substr(0, pos);
            string value = line.substr(pos + 1);
            mp[key] = value;
        }
        myfile.close();
    }
    else {
        cout << "Cannot open the file!" << endl;
    }
    return mp;

};


int main()
{
    string filename;
    filename = "config.txt";
    map<string, string> mp = read_config(filename);
    double abs_er, rel_er, x0, x1, y0, y1;
    int m, num_of_threads;
    if (mp.size() != 0) {
        abs_er = get_param<double>("absol_er", mp);
        rel_er = get_param<double>("rel_er", mp);
        x0 = get_param<double>("x0", mp);
        x1 = get_param<double>("x1", mp);
        y0 = get_param<double>("y0", mp);
        y1 = get_param<double>("y1", mp);
        m = get_param<int>("m", mp);
        num_of_threads = get_param<int>("threads", mp);
        cout << "ASD: " << num_of_threads << endl;
        thread threads[num_of_threads];
        double pr = 1E-3;

        double integral = 0;
        double interval_x = (x1 - x0) / num_of_threads;
        double x = x0;
        cout << "  Calculating...\n" << endl;
        double step1 = 1E-3;
        double step2 = step1 / 2.0;
        double integral1 = 0, integral2 = 0;
        double j = x, l = x;

        while (j < x1) {
            integral1 += integration(j, j + step1, y0, y1, m, pr);
            j += step1;
        }

        while (l < x1) {
            integral2 += integration(l, l + step2, y0, y1, m, pr);
            l += step2;
        }

        double abs_dif = abs(integral1 - integral2);
        double rel_dif = abs((integral1 - integral2) / max(integral1, integral2));

        if (abs_dif <= abs_er)
            cout << "| Absolute error is okay\t";
        else
            cout << "| Absolute error is not okay\t";

        cout << abs_dif << " vs " << abs_er << endl;

        if (rel_dif <= rel_er)
            cout << "| Relative error is okay\t";
        else
            cout << "| Relative error is not okay\t";
        cout << rel_dif << " vs " << rel_er << endl;

        auto start_time_reading = get_current_time_fenced();
        for (int i = 0; i < num_of_threads; ++i) {
            threads[i] = thread(integrateWithThreads, x, x + interval_x, y0, y1, m, pr, &integral);
            x += interval_x;
        }

        for (int i = 0; i < num_of_threads; ++i) threads[i].join();

        auto finish_time = get_current_time_fenced();
        auto total_time = finish_time - start_time_reading;

        cout << " -------------------------------------\n| Time: " << to_us(total_time)
             << " ms\n -------------------------------------" << endl;
        ofstream result, additionalInfo;
        result.open("result.txt");
        additionalInfo.open("additionalInfo.txt");
              result << integral << endl;
        additionalInfo << "| Absolute error: " << abs_dif << endl;
        additionalInfo << "| Relative error: " << rel_dif << "\n|-----------------------------" << endl;
        additionalInfo << "| Total time: " << total_time.count() << " ms\n -----------------------------" << endl;

        additionalInfo << "\t|  Threads counted: " << integral << endl;
    }
    return 0;
}
