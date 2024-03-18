在主程序的开头，除了use module之后，引入MPI库  

``` fortran
include 'mpif.h'
integer:: processor_number, my_rank, ierr
```


现在这些程序可以放在主程序变量命名后面以示MPI初始化.  

``` fortran
call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,processor_number,ierr)
```

call MPI_Comm_rank 用于获取当前进程的排名（或称为进程号）, ierr 同样是用于存储 MPI 函数的返回状态, 如果函数调用成功，ierr 的值将是 0,  processor_number 相当于获取通信域中的进程数目


在主程序的最后, deallocate的后面结束MPI进程:  
```fortran
call mpi_barrier(mpi_comm_world,ierr)
call MPI_FINALIZE(ierr)
```


------

先编译后运行:  
``` fortran
mpif90 -o xxxx.f90 xxx
```

```fortran
mpirun -np N ./xxx
```

N 就是你选择的核心数量

##### reference:

Dongxu Liu, Siyi Yang
