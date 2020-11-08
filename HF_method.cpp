#include<iostream>
#include<vector>
#include<cmath>


using namespace std;

//将特征值存放在eigvl矢量中，将特征矢量矩阵存放在eigvt中
vector<vector<double>>* eigvt = new vector<vector<double>>;
vector<double>* eigvl = new vector<double>;

//声明determinant函数（求行列式）
double determinant(vector<vector<double>>);                                 

//计算矩阵v的p行q列的元素的代数余子式，注意：此函数与determinant函数相互递归嵌套
double algCmp(vector<vector<double>> v, int p, int q)
{
    int size = v.size();
    vector<vector<double>> tmp;
    vector<double> element;
    for(int i = 0; i < size; i++)
    {
        if(i != p)
        {
            element.clear();
            for(int j = 0; j < size; j++)
            {
                if(j != q)
                {
                    element.push_back(v[i][j]);
                }
            }
            tmp.push_back(element);
        }  
    }
    return pow(-1,p+q)*determinant(tmp);
}

//计算矩阵的行列式，注意：此函数与algCmp函数相互递归嵌套
double determinant(vector<vector<double>> v)
{
    int size = v.size();
    double result = 0;
    if(size > 1)
    {
        for(int i = 0; i < size; i++)
        {
            result += v[0][i]*algCmp(v,0,i);
        }
        return result;
    }
    else
    {
        return v[0][0];
    }
}

//求矩阵的逆
vector<vector<double>>* inverse(vector<vector<double>> v)
{
    vector<vector<double>>* obj = new vector<vector<double>>;
    vector<double> element;
    double d = determinant(v);
    int size = v.size();
    for(int i = 0; i < size; i++)
    {
        element.clear();
        for(int j = 0; j < size; j++)
        {
            element.push_back(algCmp(v,j,i)/d);
        }
        obj->push_back(element);
    }
    return obj;
}

//Jacobi_method求实对称矩阵的本征值和本征矢量
void jacobi(vector<vector<double>> v)
{
	int size = v.size();
	
	vector<double> ep, eq, element, cp, cq;
	ep.resize(size, 0); eq.resize(size, 0); cp.resize(size, 0); cq.resize(size, 0);
	int p, q;
	double max = 0, err = 0.0000000000001;
	double c, s, theta;
    for(int i = 0; i < size; i++)
    {
        element.clear();
        element.resize(size,0);
        element[i] = 1;
        eigvt->push_back(element);
    }
     
	//Jacobi循环求解本征值和本征矢量
	while (true)
	{
		max = 0;
		for (int i = 0; i < size; i++)
			for (int j = i + 1; j < size; j++)
				if (abs(v[i][j]) > abs(max))
				{
					max = v[i][j];
					p = i;
					q = j;
				}
		if (abs(max) <= err)
			break;
		else
		{
			if (v[p][p] == v[q][q]) theta = atan(1);
			else theta = atan(2 * v[p][q] / (v[q][q] - v[p][p])) / 2;
			c = cos(theta);
			s = sin(theta);
			for (int i = 0; i < size; i++)
			{
				if (i == p)
				{
					ep[i] = c * c * v[p][p] - 2 * s * c * v[p][q] + s * s * v[q][q];
					eq[i] = 0;
				}
				else if (i == q)
				{
					eq[i] = s * s * v[p][p] + 2 * s * c * v[p][q] + c * c * v[q][q];
					ep[i] = 0;
				}
				else
				{
					ep[i] = c * v[p][i] - s * v[q][i];
					eq[i] = s * v[p][i] + c * v[q][i];
				}
			}
			for (int i = 0; i < size; i++)
			{
				v[p][i] = ep[i];
				v[i][p] = ep[i];
				v[q][i] = eq[i];
				v[i][q] = eq[i];
			}	

            //求出特征矢量矩阵
            for(int i = 0; i < size; i++)
            {
                cp[i] = c*(*eigvt)[i][p]-s*(*eigvt)[i][q];
                cq[i] = s*(*eigvt)[i][p]+c*(*eigvt)[i][q];
            }
            for(int i = 0; i < size; i++)
            {
                (*eigvt)[i][p] = cp[i];
                (*eigvt)[i][q] = cq[i];
            }
		}
	}

	for (int i = 0; i < size; i++)
	{
		eigvl->push_back(v[i][i]);
	}
}
 
int main()
{
    vector<vector<double>>* inv;
    
    vector<vector<double>> t;
	vector<double> e1, e2, e3, e4;
	e1.push_back(4);
	e1.push_back(-30);
	e1.push_back(60);
	e1.push_back(-35);
	e2.push_back(-30);
	e2.push_back(300);
	e2.push_back(-675);
	e2.push_back(420);
	e3.push_back(60);
	e3.push_back(-675);
	e3.push_back(1620);
	e3.push_back(-1050);
	e4.push_back(-35);
	e4.push_back(420);
	e4.push_back(-1050);
	e4.push_back(700);
	t.push_back(e1);
	t.push_back(e2);
	t.push_back(e3);
	t.push_back(e4);
	jacobi(t);
    for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << (*eigvt)[i][j] << "\t";
		}
		cout << endl;
	}
    for (int i = 0; i < 4; i++)
    cout<<(*eigvl)[i]<<endl;

    delete eigvl;
    delete eigvt;
    delete inv;
    return 0;
}