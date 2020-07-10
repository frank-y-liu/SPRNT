      subroutine tridv(node,node1,node2,coef,rank)
      real node(10),node1(10),node2(10),coef
      integer rank
      real s(3),coef1,temp
      integer t(3)
      coef1=1.0-coef
      s(1)=(node(3)-node(5))**2+(node(4)-node(6))**2
      s(2)=(node(5)-node(7))**2+(node(6)-node(8))**2
      s(3)=(node(3)-node(7))**2+(node(4)-node(8))**2
      t(1)=1
      t(2)=2
      t(3)=3
      do 10 i=1,2
        do 10 j=i+1,3
          if(s(i).lt.s(j)) then
            temp=t(i)
            t(i)=t(j)
            t(j)=temp
          end if
10    continue
      if(t(rank).eq.1)then
        node1(3)=coef*node(3)+coef1*node(5)
        node1(4)=coef*node(4)+coef1*node(6)
        node1(5)=node(5)
        node1(6)=node(6)
        node1(7)=node(7)
        node1(8)=node(8)
        node2(3)=node1(3)
        node2(4)=node1(4)
        node2(5)=node(7)
        node2(6)=node(8)
        node2(7)=node(3)
        node2(8)=node(4)
      else if(t(rank).eq.2) then
        node1(3)=coef*node(5)+coef1*node(7)
        node1(4)=coef*node(6)+coef1*node(8)
        node1(5)=node(7)
        node1(6)=node(8)
        node1(7)=node(3)
        node1(8)=node(4)
        node2(3)=node1(3)
        node2(4)=node1(4)
        node2(5)=node(3)
        node2(6)=node(4)
        node2(7)=node(5)
        node2(8)=node(6)
      else
        node1(3)=coef*node(3)+coef1*node(7)
        node1(4)=coef*node(4)+coef1*node(8)
        node1(5)=node(3)
        node1(6)=node(4)
        node1(7)=node(5)
        node1(8)=node(6)
        node2(3)=node1(3)
        node2(4)=node1(4)
        node2(5)=node(5)
        node2(6)=node(6)
        node2(7)=node(7)
        node2(8)=node(8)
      end if
      node1(9)=coef*node(9)
      node2(9)=coef1*node(9)
      end
