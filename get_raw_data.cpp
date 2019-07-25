#include"get_raw_data.h"

#define DELTA_T         0.1 //s
#define DATASET_NUM     50 //采集数据数量
#define DATA_NUM        9 //每个数据集里的数据个数
#define c cos
#define s sin
#define ITER_STEP       1e-5
#define ITER_CNT        100


/*************************************Global Parameters*************************************/
VectorXd     params_axis(4, 1);
VectorXd     params_pos(6, 1);
Vector3f     j1, j2;
Vector3f     o1, o2;
RowVector3f  g1, g2;
RowVector3f  g_dot1, g_dot2;
RowVector3f  a1, a2;
float        **imu_raw_data_1;
float        **imu_raw_data_2;
float        **imu_raw_data_online1;
float        **imu_raw_data_online2;
static float prev_angle_gyr, prev_angle_acc_gyr;

/****************************************Functions******************************************/

//函数指针声明
typedef void (*func_ptr)(const MatrixXd &input, const VectorXd &params, VectorXd &output);

float **getData( char filename[], int NUM )
{
    /*
    **  从IMU传感器获取数据到代码   
    */
    float   acc[500][3];
    float   vel[500][3];
    float   vel_dot[500][3]; 
    /*
    **  动态建立输出二维数组
    */
    float   **raw_data = (float **)malloc( NUM * sizeof(float *) );
    for( int i = 0; i < NUM; i++)
    {
        raw_data[i] = (float *)malloc( NUM * sizeof( float ) );
    }

    char    buf[300];
    char    path[] = "./";
    //char    filename[] = "rawdata.txt";
    char    path_filename[300];

    strcpy(path_filename, path);
    strcat(path_filename, filename);

    FILE *fp = fopen(path_filename, "r");

    /*
    **  定义15个文件包含的数据名
    */
    float   addr, time, tepo,
            acc_x, acc_y, acc_z, 
            vel_x, vel_y, vel_z, 
            angle_x, angle_y, angle_z;

    int     hx, hy, hz;

    if( fp )
    {
        /*
        **  不需要前两行的title
        **  然后存入buf
        */
        for( int i = 0; i < 2; i++ )
        {
            fgets( buf, sizeof(buf), fp );
            //printf("\n%s", buf);
        }
        /*
        **  读取剩下的数据
        */
        for( int i = 0; i < NUM; i++ )
        {
            fscanf(fp, "%x\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d", 
            &addr, &time, &acc_x, &acc_y, &acc_z, &vel_x, &vel_y, &vel_z,
            &angle_x, &angle_y, &angle_z, &tepo, &hx, &hy, &hz );
            //printf("\n%x\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%d\n%d\n%d\n",addr, time, acc_x, acc_y, acc_z, vel_x, vel_y, vel_z,
            //angle_x, angle_y, angle_z, tepo, hx, hy, hz);

            acc[i][0] = acc_x;
            acc[i][1] = acc_y;
            acc[i][2] = acc_z;
            
            vel[i][0] = vel_x;
            vel[i][1] = vel_y;
            vel[i][2] = vel_z;
            
        }
        fclose(fp);
        /*
        **  计算角速度在x,y,z的微分
        */
        for( int i = 0; i < NUM; i++ )
        {
            for( int j = 0; j < 3; j++ )
            {
                if( i > 1 )
                {
                    vel_dot[i][j] = ( vel[i-2][j] - 8 * vel[i-1][j] + 8 * vel[i+1][j] - vel[i+2][j] ) / 12 * DELTA_T;
                }
                else
                {
                    vel_dot[i][j] = ( 8 * vel[i+1][j] - vel[i+2][j] ) / 12 * DELTA_T;
                }
            }
            
            /*
            **  把三个数据融合到一个数据输出
            */
            int k = 0;
            for( int j = 0; j < 3; j++ )
            {
                raw_data[i][k] = acc[i][j];
                raw_data[i][k+3] = vel[i][j]; 
                raw_data[i][k+6] = vel_dot[i][j];
                k++;
                 
            }
            
            /*
            **  打印输出数组
            */
            // for( int k = 0; k < DATA_NUM; k++ )
            // {
            //     printf("%f\n",raw_data[i][k] );
            // }  
              
        }  
        return raw_data;
    }
    
    else
    {
        /*
        **  如果是空指针，程序停止
        */
        printf("cannot open %s", filename);
        return NULL;
    }

}

void get_raw_data()
{
    /*
    **  获取两个imu的数据,这里只需要角速度
    */
    char imu_filename_1[] = "rawdata1.txt";
    char imu_filename_2[] = "rawdata2.txt";

    imu_raw_data_1 = getData( imu_filename_1, DATASET_NUM );
    imu_raw_data_2 = getData( imu_filename_2, DATASET_NUM );
    
    /*
    **  prinf
    */
    // for( int i = 0; i < 50; i++ )
    // {
    //     for( int j = 0; j < 9; j++ )
    //     {
    //         printf("%f ", imu_raw_data_1[i][j] );
    //     }
    //     printf("\n");
    // }
    
}

void get_pos(const MatrixXd &input, const VectorXd &params, VectorXd &output)
{
    /*
    **  获取关节相对于两个imu的位置
    */

    //定义6个待求参数
    float o1x = params(0, 0);
    float o1y = params(1, 0);
    float o1z = params(2, 0);
    float o2x = params(3, 0);
    float o2y = params(4, 0);
    float o2z = params(5, 0);
  
    for( int i = 0; i < input.rows(); i++ )
    {

        //角加速度计算值
        float acc_joint1_x = input(i, 4) * ( input(i, 3) * o1y - input(i, 4) * o1x ) - input(i, 5) * ( input(i, 5) * o1x - input(i, 3) * o1z ) + ( input(i, 7) * o1z - input(i, 8) * o1y );
        float acc_joint1_y = input(i, 5) * ( input(i, 4) * o1z - input(i, 5) * o1y ) - input(i, 3) * ( input(i, 3) * o1y - input(i, 4) * o1x ) + ( input(i, 8) * o1x - input(i, 6) * o1z );
        float acc_joint1_z = input(i, 3) * ( input(i, 5) * o1x - input(i, 3) * o1z ) - input(i, 4) * ( input(i, 4) * o1z - input(i, 5) * o1y ) + ( input(i, 6) * o1y - input(i, 7) * o1x );

        float acc_joint2_x = input(i, 13) * ( input(i, 12) * o2y - input(i, 13) * o2x ) - input(i, 14) * ( input(i, 14) * o2x - input(i, 12) * o2z ) + ( input(i, 16) * o2z - input(i, 17) * o2y );
        float acc_joint2_y = input(i, 14) * ( input(i, 13) * o2z - input(i, 14) * o2y ) - input(i, 12) * ( input(i, 12) * o2y - input(i, 13) * o2x ) + ( input(i, 17) * o2x - input(i, 15) * o2z );
        float acc_joint2_z = input(i, 12) * ( input(i, 14) * o2x - input(i, 12) * o2z ) - input(i, 13) * ( input(i, 13) * o2z - input(i, 14) * o2y ) + ( input(i, 15) * o2y - input(i, 16) * o2x );

        //目标函数
        output(i, 0) =  sqrt(
                            pow( (input(i, 0) - acc_joint1_x ), 2) + pow( (input(i, 1) - acc_joint1_y ), 2) + pow( (input(i, 2) - acc_joint1_z ), 2)
                            ) - 
                        sqrt(
                            pow( (input(i, 9) - acc_joint2_x ), 2) + pow( (input(i, 10) - acc_joint2_y ), 2) + pow( (input(i, 11) - acc_joint2_z ), 2)
                        );
    }

}

void get_axis(const MatrixXd &input, const VectorXd &params, VectorXd &output)
{
    /*
    **  获取关节相对于两个imu的方向
    */

    float theta_1 = params(0, 0);
    float theta_2 = params(1, 0);
    float phi_1   = params(2, 0);
    float phi_2   = params(3, 0);

    for( int i = 0; i < input.rows(); i++ )
    {
        //目标模型
        output(i, 0) =  sqrt(
                            pow( (input(i, 1) * s(theta_1) - input(i, 2) * c(phi_1) * s(theta_1)), 2) +
                            pow( (input(i, 2) * c(phi_1) * c(theta_1) - input(i, 0) * s(phi_1)), 2) +
                            pow( (input(i, 0) * c(phi_1) * s(theta_1) - input(i, 1) * c(phi_1) * c(theta_1)), 2) ) -
                        sqrt(
                            pow( (input(i, 4) * s(theta_2) - input(i, 5) * c(phi_2) * s(theta_2)), 2) +
                            pow( (input(i, 5) * c(phi_2) * c(theta_2) - input(i, 3) * s(phi_2)), 2) +
                            pow( (input(i, 3) * c(phi_2) * s(theta_2) - input(i, 4) * c(phi_2) * c(theta_2)), 2) );
    }
}

void get_jacobian(  func_ptr func,
                    const MatrixXd &input,
                    const VectorXd &params,
                    MatrixXd &output
                    )
{
    /*
    **  获取高斯牛顿法迭代式子里的Jacobian
    */
    int m = input.rows();   //数据数量
    int n = params.rows();  //未知参数数量

    VectorXd out0(m, 1);
    VectorXd out1(m, 1);
    VectorXd param0(n, 1);
    VectorXd param1(n, 1);

    for(int j = 0; j < n; j++ )
    {
        param0 = params;
        param1 = params;
        // cout << param0 << "\n" << endl;
        param0(j, 0) -= ITER_STEP;
        param1(j, 0) += ITER_STEP;
        func(input, param0, out0);
        func(input, param1, out1);
       
        output.block(0, j, m, 1) = (out1 - out0) / ( 2 * ITER_STEP );
        // cout << output << "\n" << endl;
    }
}

void gauss_newton(  func_ptr func,
                    const MatrixXd &input,
                    const VectorXd &output,
                    VectorXd &params
                    )
{
    /*
    **  高斯牛顿法
    */
    int m = input.rows();
    int n = params.rows();

    //jacobian
    MatrixXd jmat(m, n);
    VectorXd r(m, 1);
    VectorXd tmp(m, 1);

    float pre_mse = 0.0;
    float mse;
    
    for(int i = 0; i < ITER_CNT; i++ )
    {
        mse = 0.0;
        func( input, params, tmp );
        r = output - tmp;
        get_jacobian(func, input, params, jmat);
    
        //均方误差
        mse = r.transpose() * r;
        mse /= m;
        if( fabs(mse - pre_mse) < 1e-8)
        {
            break;
        }
        pre_mse = mse;

        //参数更新
        VectorXd delta = (jmat.transpose() * jmat).inverse() * jmat.transpose() * r;
        printf("i = %d, mse %lf \n", i, mse);
        params += delta;
        
    }

    cout << "params:" << params.transpose() << endl;
}   

int imu_joint_pos_data_fit()
{
    /*
    **  计算关节位置的数据输入接口
    */

    MatrixXd input(DATASET_NUM, 18);
    VectorXd output(DATASET_NUM, 1);
    
    for( int i = 0; i < DATASET_NUM; i++ )
    {
        int k = 0;
        for( int j = 0; j < 9; j++ )
        {   
            input(i, k) = imu_raw_data_1[i][j];
            k++;
        }
        for( int j = 0; j < 9; j++ )
        {
            input(i, k) = imu_raw_data_2[i][j];
            k++;
        }
        output(i,0) = 0;
    }
    
    params_pos << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    gauss_newton(get_pos, input, output, params_pos);

    o1 << params_pos(0, 0), params_pos(1, 0), params_pos(2, 0);
    o2 << params_pos(3, 0), params_pos(4, 0), params_pos(5, 0);

    return 0;
    
}

int imu_joint_axis_data_fit()
{
    /*
    **  计算关节轴向的数据输入接口
    */

    MatrixXd input(DATASET_NUM, 6);
    VectorXd output(DATASET_NUM, 1);
    
    for( int i = 0; i < DATASET_NUM; i++ )
    {
        int k = 0;
        for( int j = 3; j < 6; j++ )
        {   
            input(i, k) = imu_raw_data_1[i][j];
            k++;
        }
        for( int j = 3; j < 6; j++ )
        {
            input(i, k) = imu_raw_data_2[i][j];
            k++;
        }
        output(i,0) = 0;
    }
    
    params_axis << 0.5, 0.5, 0.5, 0.5;
    
    gauss_newton(get_axis, input, output, params_axis);
    
    j1 << c(params_axis(2, 0)) * c(params_axis(0, 0)), c(params_axis(2, 0)) * s(params_axis(0, 0)), s(params_axis(2, 0));
    j2 << c(params_axis(3, 0)) * c(params_axis(1, 0)), c(params_axis(3, 0)) * s(params_axis(1, 0)), s(params_axis(3, 0));

    return 0;
    
}

float get_angle_acc(Vector3f j1, Vector3f j2,
                    RowVector3f a1, RowVector3f a2,
                    RowVector3f g1, RowVector3f g2,
                    RowVector3f g_dot1, RowVector3f g_dot2,
                    Vector3f o1, Vector3f o2
                    )
{
    /*
    **  计算基于imu加速度数据解得的角度
    */
    Vector3f        c, x1, x2, y1, y2;
    Vector2f        acc1, acc2;
    RowVector3f     a1_dot, a2_dot;

    float p1, p2, q1, q2;
    float angle_acc;

    c << 1, 0, 0;

    /*
    **  处理o1, o2
    */
    o1 = o1 - j1 * ( o1.dot(j1) + o2.dot(j2) ) / 2;
    o2 = o2 - j2 * ( o1.dot(j1) + o2.dot(j2) ) / 2;

    a1_dot = a1 - ( g1.cross( g1.cross(o1) ) + g_dot1.cross(o1) );
    a2_dot = a2 - ( g2.cross( g2.cross(o2) ) + g_dot2.cross(o2) );

    x1 = j1.cross(c);
    y1 = j1.cross(x1);
    x2 = j2.cross(c);
    y2 = j2.cross(x2);

    p1 = a1_dot * x1;
    p2 = a1_dot * y1;
    q1 = a2_dot * x2;
    q2 = a2_dot * y2;

    acc1 << p1, p2;
    acc2 << q1, q2;

    angle_acc = acos( acc1.dot(acc2) / (acc1.norm() * acc2.norm()) );

    return angle_acc;
}

void test_angle()
{
    float angle_acc, angle_gyr, angle_acc_gyr;
    float sum = 0;         
    int   cnt = 0;
    float lambda = 0.01;

    // for test
    char imu_dataonline1[] = "rawdata_online1.txt";
    char imu_dataonline2[] = "rawdata_online2.txt";
    imu_raw_data_online1 = getData( imu_dataonline1, 500 );
    imu_raw_data_online2 = getData( imu_dataonline2, 500 );
    
    FILE *fp;
    fp = fopen("data.txt","w");
    if(fp == NULL)
    {
        printf("File cannot open");
        exit(0);
    }


    for( int i = 0; i < 500; i++ )
    {
        cnt++;

        a1 << imu_raw_data_online1[i][0], imu_raw_data_online1[i][1], imu_raw_data_online1[i][2];
        a2 << imu_raw_data_online2[i][0], imu_raw_data_online2[i][1], imu_raw_data_online2[i][2];
        g1 << imu_raw_data_online1[i][3], imu_raw_data_online1[i][4], imu_raw_data_online1[i][5];
        g2 << imu_raw_data_online2[i][3], imu_raw_data_online2[i][4], imu_raw_data_online2[i][5];
        g_dot1 << imu_raw_data_online1[i][6], imu_raw_data_online1[i][7], imu_raw_data_online1[i][8];
        g_dot2 << imu_raw_data_online2[i][6], imu_raw_data_online2[i][7], imu_raw_data_online2[i][8];

        angle_acc = get_angle_acc( j1, j2, a1, a2, g1, g2, g_dot1, g_dot2, o1, o2 );

        sum = sum + g1 * j1 - g2 * j2;

        if( cnt > 3 )
        {
            /*
            **  计算基于imu陀螺仪数据所解得的角度
            */ 
            angle_gyr = sum * DELTA_T;

            /*
            **  互补算法融合两种不同方法算出来的角度
            */
            angle_acc_gyr = lambda * angle_acc + (1-lambda) * ( prev_angle_acc_gyr + angle_gyr - prev_angle_gyr );
            cout << "angle: " << angle_acc_gyr << endl;
            cnt = 0;
            /*
            **  数据写入文档
            */
            fprintf(fp,"%f\n", angle_acc_gyr);
        }
        /*
        **  数据更新
        */
        prev_angle_acc_gyr = angle_acc_gyr;
        prev_angle_gyr = angle_gyr;

    }
    fclose(fp);
}

int main()
{
    get_raw_data();
    imu_joint_axis_data_fit();
    imu_joint_pos_data_fit();

    /*
    **  refers to an arbitrary point along the joint axis
    **  we shift it as close as possible to the sensors by applying:
    */
    o1 = o1 - j1 * ( o1.dot(j1) + o2.dot(j2) ) / 2;
    o2 = o2 - j2 * ( o1.dot(j1) + o2.dot(j2) ) / 2;
    
    test_angle();
    
    return 0;
    
}
