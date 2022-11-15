1.使用libigl将Obj先转化为off格式  √
2.使用CGAL的package直接将off读入到ployhedron结构中
3.使用ployhedron进行半边操作
4.将ployhedron的数据结构转化为eigen格式
5.用libigl显示 
	根据libigl example的模板，fetchcontent
	将fetchcontent_source_dir_libigl="path-to-libigl"
	然后generate就行
	将Igl从include目录中移除，要不会包含它而不是下载的子模块

6. GetMeanCuravaterOpreator()问题
	1.He为何是向外的法向，即给更新顶点会增加噪声
	2.He过大，超出屏幕