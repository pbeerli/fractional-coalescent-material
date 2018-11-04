
void int_to_binary(int x, char* val)
{
  long cnt, mask = 1 << 31;
  long z=0;
  //  if(x>mask)
  //{
  //  error("int_to_binary() failed because the parameter is larger than the mask");
  //}
  for(cnt=1;cnt<=32;++cnt)
    {
      val[z++]=(((x & mask) == 0) ? '0' : '1');
      x <<= 1;
    }
}
