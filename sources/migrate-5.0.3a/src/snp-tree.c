/*

SNP tree changing algorithms: this assumes that we will have maximally two alleles per SNP in a population.
It further assumes that the SNP are unlinked or that we have a short sequence around linked SNPs. 

(c) Beerli 2006, Tallahassee

*/

  void change_snp_tree(world_fmt *world);




///
/// changes a SNP tree. Such a tree has maximally alleles * populations tips, where each tree has also
/// additional events that show the MRCA for each population.
/// A-O--------* 
///            |
///            *---
///            |
/// C----O-----*
///
/// A1-O----*--* 
///         |  |
/// A2--O-|-*  *---
///            |
/// C1---O-----*
///

  void change_snp_tree(world_fmt *world)
{

}





