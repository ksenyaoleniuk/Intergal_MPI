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

//Моя формула закоментована, оскільки ми робили з Оксаною, ми використовували її функцію.
double func_calculation(int m, double x1, double x2) {
//    double sum1 = 0;
//    double sum2 = 0;
//    double g;
//    for (int i = 1; i <= m; ++i)
//    {
//        sum1 = c(i)*(exp(-1/pi·sum((x-A(i))^2))·cos(pi·sum((x-A(i))^2)))
//    }
//
//    g = - sum1 * sum2;
//
//    return g;
    double g, sum;
    int j;
    for (int i = -2; i <= 2; ++i) {
        j = i;
        sum += 1 / (5 * (i + 2) + j + 3 + pow(x1 - 16 * j, 6) + pow(x2 - 16 * i, 6));
    }
    g = pow(0.002 + sum, -1);
    return g;
}

double integrate(double x0, double x, double y0, double y, int m, double pr) {
    assert (m >= 5);
    double sum = 0;
    for (double i = x0; i <= x; i += pr) {
        for (double j = y0; j <= y; j += pr) {
            sum += func_calculation(m, i + pr / 2.0, j + pr / 2.0) * pr * pr;
        }
    }
    return sum;
}


template<class T>
T get_param(string key, map<string, string> myMap) {
    istringstream ss(myMap[key]);
    T val;
    ss >> val;
    return val;
}

void thread_integration(double x0, double x, double y0, double y, int m, double pr, double *r) {
    auto answer = integrate(x0, x, y0, y, m, pr);
    lock_guard<mutex> lg(mx);
    *r += answer;
}

map<string, string> read_config(string filename) {
    string line, delimiter;
    ifstream myfile(filename);
    map<string, string> mp;

    delimiter = "=";


    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            int pos = line.find(delimiter);
            string key = line.substr(0, pos);
            string value = line.substr(pos + 1);
            mp[key] = value;
        }

        myfile.close();
    } else {
        cout << "Error with opening the file!" << endl;
    }
    return mp;

};

template<class D>
inline long long to_us(const D &d) {
    return std::chrono::duration_cast<chrono::microseconds>(d).count();
}


inline std::chrono::high_resolution_clock::time_point get_current_time_fenced() {
    std::atomic_thread_fence(memory_order_seq_cst);
    auto res_time = std::chrono::high_resolution_clock::now();
    std::atomic_thread_fence(memory_order_seq_cst);
    return res_time;
}


int main() {
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Get_processor_name(procname, &len);

    string filename;
//    cout << "Please enter name of configuration file with extension '.txt':";
//    cin >> filename;

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

        thread threads[num_of_threads];
        double pr = 1E-3;

        double integral = 0;
        double interval_x = (x1 - x0) / num_of_threads;
        double x = x0;
        cout << "  \t\t\t\t\t\t\t\t\t\t\t\t\tGetting result\n" << endl;
        double step1 = 1E-3;
        double step2 = step1 / 2.0;
        double firtInt = 0, secondInt = 0;
        double j = x, l = x;

        while (j < x1) {
            firtInt += integrate(j, j + step1, y0, y1, m, pr);
            j += step1;
        }

        while (l < x1) {
            secondInt += integrate(l, l + step2, y0, y1, m, pr);
            l += step2;
        }

        double abs_dif = abs(firtInt - secondInt);
        double rel_dif = abs((firtInt - secondInt) / max(firtInt, secondInt));


        auto start_time_reading = get_current_time_fenced();
        for (int i = 0; i < num_of_threads; ++i) {
            threads[i] = thread(thread_integration, x, x + interval_x, y0, y1, m, pr, &integral);
            x += interval_x;
        }

        for (int count = 0; count < num_of_threads; ++count) threads[count].join();

        auto ending_time = get_current_time_fenced();
        auto total_time = ending_time - start_time_reading;

        cout << "Absolute error: real vs desirable - " << abs_dif << " vs " << abs_er << endl;
        cout << "Relative error: real vs desirable - " << rel_dif << " vs " << rel_er << endl;
        cout << "Result: " << integral << endl;
        cout << "Total time: " << to_us(total_time) << endl;

        ofstream result, additionalInf;
        result.open("result.txt");
        additionalInf.open("additionInf.txt");
        additionalInf << "Integral: " << integral << endl;
        additionalInf << "Absolute error: " << abs_dif << endl;
        additionalInf << "Relative error: " << rel_dif << endl;
        result << total_time.count() << endl;
    }
    return 0;
}