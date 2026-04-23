#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"

// 如果你不需要使用 matrix 类，请将 IGNORE_MATRIX 改为 0
// #define IGNORE_MATRIX 0
#define IGNORE_MATRIX 1

#if IGNORE_MATRIX

class matrix {
private:

    // m行n列的矩阵，用动态二维数组存储，每个元素是分数类实例
    int m, n;
    fraction **data;

public:

    // 默认构造函数
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    // 构造函数，构建 m_*n_ 的矩阵，矩阵元素设为0。
    matrix(int m_, int n_) : m(m_), n(n_) {
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = fraction(0);
            }
        }
    }

    // 拷贝构造函数，构建与 obj 完全相同的矩阵。
    matrix(const matrix &obj) : m(obj.m), n(obj.n) {
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = obj.data[i][j];
            }
        }
    }

    // 移动拷贝构造函数。
    matrix(matrix &&obj) noexcept : m(obj.m), n(obj.n), data(obj.data) {
        obj.m = obj.n = 0;
        obj.data = nullptr;
    }

    // 析构函数。
    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    // 重载赋值号。
    matrix &operator=(const matrix &obj) {
        if (this != &obj) {
            if (data != nullptr) {
                for (int i = 0; i < m; i++) {
                    delete[] data[i];
                }
                delete[] data;
            }
            
            m = obj.m;
            n = obj.n;
            data = new fraction*[m];
            for (int i = 0; i < m; i++) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; j++) {
                    data[i][j] = obj.data[i][j];
                }
            }
        }
        return *this;
    }

    // 重载括号，返回矩阵的第i行(1-based)、第j列(1-based)的元素的引用。如果 i、j 不合法，抛出 matrix_error 错误。
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 1 || j > n) {
            throw matrix_error();
        }
        return data[i-1][j-1];
    }
    
    // const版本的重载括号
    const fraction &operator()(int i, int j) const {
        if (i < 1 || i > m || j < 1 || j > n) {
            throw matrix_error();
        }
        return data[i-1][j-1];
    }

    // 重载乘号，返回矩阵乘法 lhs * rhs 的结果。如果 lhs 的列数与 rhs 的行数不相等，抛出 matrix_error 错误。
    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m) {
            throw matrix_error();
        }
        
        matrix result(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; i++) {
            for (int j = 0; j < rhs.n; j++) {
                fraction sum(0);
                for (int k = 0; k < lhs.n; k++) {
                    sum = sum + lhs.data[i][k] * rhs.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }
        return result;
    }

    // 返回矩阵的转置。若矩阵为空，抛出 matrix_error 错误。
    matrix transposition() {
        if (m == 0 || n == 0) {
            throw matrix_error();
        }
        
        matrix result(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    // 返回矩阵的行列式。建议用高斯消元实现。若矩阵不是方阵或为空，抛出 matrix_error 错误。
    fraction determination() {
        if (m != n || m == 0) {
            throw matrix_error();
        }
        
        matrix temp(*this);
        fraction det(1);
        
        for (int i = 0; i < n; i++) {
            // 找到主元
            int pivot = i;
            for (int j = i + 1; j < n; j++) {
                if (!(temp(j + 1, i + 1) == fraction(0)) && 
                    (temp(pivot + 1, i + 1) == fraction(0))) {
                    pivot = j;
                }
            }
            
            if (temp(pivot + 1, i + 1) == fraction(0)) {
                return fraction(0);
            }
            
            // 交换行
            if (pivot != i) {
                for (int j = 0; j < n; j++) {
                    fraction tmp = temp.data[i][j];
                    temp.data[i][j] = temp.data[pivot][j];
                    temp.data[pivot][j] = tmp;
                }
                det = fraction(0) - det; // 改变符号
            }
            
            det = det * temp(i + 1, i + 1);
            
            // 消元
            for (int j = i + 1; j < n; j++) {
                if (!(temp(j + 1, i + 1) == fraction(0))) {
                    fraction factor = temp(j + 1, i + 1) / temp(i + 1, i + 1);
                    for (int k = i; k < n; k++) {
                        temp.data[j][k] = temp.data[j][k] - factor * temp.data[i][k];
                    }
                }
            }
        }
        
        return det;
    }
    
    int get_rows() const { return m; }
    int get_cols() const { return n; }
};

#endif

class resistive_network {
private:
	
	
    /* 以下是建议的类成员框架，你可以选择使用或者自己实现。
    
    // 节点数量 和 接线数量
    int interface_size, connection_size;

    // 邻接矩阵A，电导矩阵C，Laplace矩阵(A^tCA) (具体定义见 reading_materials.pdf)
    matrix adjacency, conduction, laplace;
    
    */ 

    // 节点数量 和 接线数量
    int interface_size, connection_size;

    // 邻接矩阵A，电导矩阵C，Laplace矩阵(A^tCA) (具体定义见 reading_materials.pdf)
    matrix adjacency, conduction, laplace;

    // 辅助函数：求解线性方程组 Ax = b
    matrix solve_linear_system(const matrix& A, const matrix& b);

public:

    //****************************
	// 注意，请勿私自修改以下函数接口的声明！
    //****************************	

    // TODO: 设置电阻网络。节点数量为interface_size_，接线数量为connection_size_。
    //       对于 1<=i<=connection_size_，从节点from[i-1]到节点to[i-1]有接线，对应电阻为resistance[i-1]。
    //       保证接线使得电阻网络联通，from[i-1] < to[i-1]，resitance[i-1] > 0，均合法。
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]);

    ~resistive_network() = default;

    // TODO: 返回节点 interface_id1 和 interface_id2 (1-based)之间的等效电阻。
    //       保证 interface_id1 <= interface_id2 均合法。
    fraction get_equivalent_resistance(int interface_id1, int interface_id2);

    // TODO: 在给定节点电流I的前提下，返回节点id(1-based)的电压。认为节点interface_size(1-based)的电压为0。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电流为 current[i-1]。
    //       保证 current 使得电阻网络有解，id < interface_size 合法。
    fraction get_voltage(int id, fraction current[]);


    // TODO: 在给定节点电压U的前提下，返回电阻网络的功率。
    //       对于 1<=i<=interface_size，节点i (1-based) 对应电压为 voltage[i-1]。
    //       保证 voltage 合法。
    fraction get_power(fraction voltage[]);
};

// 实现构造函数
resistive_network::resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
    : interface_size(interface_size_), connection_size(connection_size_),
      adjacency(interface_size_, connection_size_),
      conduction(connection_size_, connection_size_),
      laplace(interface_size_, interface_size_) {
    
    // 构建邻接矩阵A
    for (int i = 0; i < connection_size_; i++) {
        adjacency(from[i], i + 1) = fraction(1);
        adjacency(to[i], i + 1) = fraction(-1);
    }
    
    // 构建电导矩阵C（对角矩阵）
    for (int i = 0; i < connection_size_; i++) {
        conduction(i + 1, i + 1) = fraction(1) / resistance[i];
    }
    
    // 计算Laplace矩阵 L = A^t * C * A
    matrix At = adjacency.transposition();
    matrix temp = At * conduction;
    laplace = temp * adjacency;
}

// 辅助函数：求解线性方程组 Ax = b（使用高斯消元法）
matrix resistive_network::solve_linear_system(const matrix& A, const matrix& b) {
    int n = A.get_rows();
    matrix augmented(n, n + 1);
    
    // 构建增广矩阵
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented(i + 1, j + 1) = A(i + 1, j + 1);
        }
        augmented(i + 1, n + 1) = b(i + 1, 1);
    }
    
    // 高斯消元
    for (int i = 0; i < n; i++) {
        // 找到主元
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (!(augmented(j + 1, i + 1) == fraction(0)) && 
                (augmented(pivot + 1, i + 1) == fraction(0))) {
                pivot = j;
            }
        }
        
        if (augmented(pivot + 1, i + 1) == fraction(0)) {
            throw resistive_network_error();
        }
        
        // 交换行
        if (pivot != i) {
            for (int j = i + 1; j <= n + 1; j++) {
                fraction tmp = augmented(i + 1, j);
                augmented(i + 1, j) = augmented(pivot + 1, j);
                augmented(pivot + 1, j) = tmp;
            }
        }
        
        // 归一化主元行
        fraction pivot_val = augmented(i + 1, i + 1);
        for (int j = i + 1; j <= n + 1; j++) {
            augmented(i + 1, j) = augmented(i + 1, j) / pivot_val;
        }
        
        // 消元
        for (int j = 0; j < n; j++) {
            if (j != i && !(augmented(j + 1, i + 1) == fraction(0))) {
                fraction factor = augmented(j + 1, i + 1);
                for (int k = i + 1; k <= n + 1; k++) {
                    augmented(j + 1, k) = augmented(j + 1, k) - factor * augmented(i + 1, k);
                }
            }
        }
    }
    
    // 提取解
    matrix solution(n, 1);
    for (int i = 0; i < n; i++) {
        solution(i + 1, 1) = augmented(i + 1, n + 1);
    }
    
    return solution;
}

// 计算等效电阻
fraction resistive_network::get_equivalent_resistance(int interface_id1, int interface_id2) {
    if (interface_id1 == interface_id2) {
        return fraction(0);
    }
    
    // 创建简化的Laplace矩阵（去掉interface_id2行和列）
    int n = interface_size - 1;
    matrix reduced_laplace(n, n);
    
    int row = 0, col = 0;
    for (int i = 1; i <= interface_size; i++) {
        if (i == interface_id2) continue;
        col = 0;
        for (int j = 1; j <= interface_size; j++) {
            if (j == interface_id2) continue;
            reduced_laplace(row + 1, col + 1) = laplace(i, j);
            col++;
        }
        row++;
    }
    
    // 计算行列式
    fraction det = reduced_laplace.determination();
    
    // 计算余子式
    matrix cofactor(n - 1, n - 1);
    row = col = 0;
    for (int i = 1; i <= n; i++) {
        if (i == interface_id1) continue;
        col = 0;
        for (int j = 1; j <= n; j++) {
            if (j == interface_id1) continue;
            cofactor(row + 1, col + 1) = reduced_laplace(i, j);
            col++;
        }
        row++;
    }
    
    fraction cofactor_det = cofactor.determination();
    
    return cofactor_det / det;
}

// 计算电压
fraction resistive_network::get_voltage(int id, fraction current[]) {
    // 创建简化的Laplace矩阵（去掉最后一行和列）
    int n = interface_size - 1;
    matrix reduced_laplace(n, n);
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            reduced_laplace(i, j) = laplace(i, j);
        }
    }
    
    // 创建电流向量（去掉最后一个元素，但需要调整）
    matrix current_vector(n, 1);
    for (int i = 0; i < n; i++) {
        current_vector(i + 1, 1) = current[i];
    }
    
    // 求解线性方程组
    matrix solution = solve_linear_system(reduced_laplace, current_vector);
    
    return solution(id, 1);
}

// 计算功率
fraction resistive_network::get_power(fraction voltage[]) {
    fraction total_power(0);
    
    // 计算每条支路上的功率
    for (int i = 0; i < connection_size; i++) {
        // 找到这条支路连接的两个节点
        int node1 = -1, node2 = -1;
        for (int j = 1; j <= interface_size; j++) {
            if (adjacency(j, i + 1) == fraction(1)) {
                node1 = j;
            } else if (adjacency(j, i + 1) == fraction(-1)) {
                node2 = j;
            }
        }
        
        if (node1 != -1 && node2 != -1) {
            // 计算支路电压差
            fraction voltage_diff = voltage[node1 - 1] - voltage[node2 - 1];
            
            // 计算支路电流
            fraction conductance = conduction(i + 1, i + 1);
            fraction current = conductance * voltage_diff;
            
            // 计算功率 P = I^2 * R = V^2 / R = I * V
            fraction power = current * voltage_diff;
            total_power = total_power + power;
        }
    }
    
    return total_power;
}

#endif //SRC_HPP