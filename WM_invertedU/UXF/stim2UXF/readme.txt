0. 独立运行 single 中的文件C语言工程15次，命名为dopamine_stim2_uxf* 其中*为实验次数, 都拷贝到data 下， 用于画带方差的 uxf图；
1.数据生成完成后，用UXF_2017_8_10_sparse.m进行画图

2.放电率是统计量与u，x连续量的计算不同；为了保持一致都取bin2=10ms；firerate 在一个Tbin中只有一个值，而u，x是Tbin内所有ux的均值作为一个值。