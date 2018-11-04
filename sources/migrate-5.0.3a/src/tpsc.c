// TPSC -- transition probability structured coalescence
// 
//
void calculate_pop_prob(node *p, world_fmt *world)
{
  node *left = crawlback(p->next);
  node *right = crawlback(p->next->next);

  MYREAL * leftval = left->probs;
  MYREAL * rightval = right->probs;
  
