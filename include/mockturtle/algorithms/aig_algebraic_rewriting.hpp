/*!
  \file aig_algebraic_rewriting.hpp
  \brief AIG algebraric rewriting

  EPFL CS-472 2021 Final Project Option 1
*/

#pragma once

#include "../networks/aig.hpp"
#include "../views/depth_view.hpp"
#include "../views/topo_view.hpp"
#include <kitty/print.hpp>

namespace mockturtle
{

namespace detail
{

template<class Ntk>
class aig_algebraic_rewriting_impl
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  aig_algebraic_rewriting_impl( Ntk& ntk )
    : ntk( ntk )
  {
    static_assert( has_level_v<Ntk>, "Ntk does not implement depth interface." );
  }

  void run()
  {
    bool cont{true}; /* continue trying */
    while ( cont )
    {
      cont = false; /* break the loop if no updates can be made */
      ntk.foreach_gate( [&]( node n ){
        if ( try_algebraic_rules( n ) )
        {
          ntk.update_levels();
          cont = true;
        }
      });
    }
  }

private:
  /* Try various algebraic rules on node n. Return true if the network is updated. */
  bool try_algebraic_rules( node n )
  {
    if ( try_associativity( n ) )
      return true;
    if ( try_distributivity( n ) )
      return true;
    if ( try_3_layer_distributivity( n ) )
	  return true;
    /* TODO: add more rules here... */

    return false;
  }

  /* Try the associativity rule on node n. Return true if the network is updated. */
  bool try_associativity( node n )
  { 	
  	std::vector<signal> sigVec;
  	
	ntk.foreach_fanin( n, [&]( signal sig ){ sigVec.push_back(sig); });	//Extracts the fanin signals
	
	if(sigVec.size() != 2){
		return false;	//No optimization if we don't have 2 fanin signals
	}
	
	int32_t diffInLevel = ntk.level(ntk.get_node(sigVec.at(0))) - ntk.level(ntk.get_node(sigVec.at(1)));	//Signed difference in level between fanins
	signal sigH, a;
	
	if(diffInLevel >= 2){	//We could benefit from inverting the associativity's parenthesis (sigVec[0] higher)
		sigH = sigVec.at(0);
		a = sigVec.at(1);
	}else if(diffInLevel <= -2){	//We could benefit from inverting the associativity's parenthesis (sigVec[1] higher)
		sigH = sigVec.at(1);
		a = sigVec.at(0);
	}else{
		return false;	//No gain from associativity
	}

	if(ntk.is_complemented(sigH) == true || ntk.fanout_size(ntk.get_node(sigH)) != 1){
		return false;	//Associativity requires uncomplemented node with no fanout here
	}

	//	Fanins of fanin
	sigVec.clear();
	ntk.foreach_fanin( ntk.get_node(sigH), [&]( signal sig ){ sigVec.push_back(sig); });	//Extracts the fanin signals

	if(sigVec.size() != 2){
		return false;	//No optimization if we don't have 2 fanin signals
	}

	signal b, c;
	if( ntk.level(ntk.get_node(sigVec.at(0))) >= ntk.level(ntk.get_node(sigVec.at(1))) ){	//Sets as b the higher node, to invert with a
		b = sigVec.at(0);
		c = sigVec.at(1);
	}else{
		b = sigVec.at(1);
		c = sigVec.at(0);
	}

 	bool cpla = ntk.is_complemented(a), cplb = ntk.is_complemented(b);
	if(cplb){a = !a;}
	if(cpla){b = !b;}
	
	ntk.replace_in_node(n, ntk.get_node(a), b);	//Exchanges a and b
	ntk.replace_in_node(ntk.get_node(sigH), ntk.get_node(b), a);

	return true;
  }
  
  /* Try the distributivity rule on node n. Return true if the network is updated. */
  bool try_distributivity( node n )
  {
  	std::vector<signal> sigVec;
  	
	ntk.foreach_fanin( n, [&]( signal sig ){ sigVec.push_back(sig); });	//Extracts the fanin signals
	
	if(sigVec.size() != 2){
		return false;	//No optimization if we don't have 2 fanin signals
	}
	
	if(ntk.fanout_size(ntk.get_node(sigVec.at(0))) != 1 || ntk.fanout_size(ntk.get_node(sigVec.at(1))) != 1){
		return false;	//Fanouts in middle nodes, can't modify them
	}
	
	std::vector<signal> sigVec1, sigVec2;
	ntk.foreach_fanin( ntk.get_node(sigVec.at(0)), [&]( signal sig ){ sigVec1.push_back(sig); });	//Extracts the fanin signals from first fanin
	ntk.foreach_fanin( ntk.get_node(sigVec.at(1)), [&]( signal sig ){ sigVec2.push_back(sig); });	//Extracts the fanin signals from second fanin
	
	if(sigVec1.size() != 2 || sigVec2.size() != 2){
		return false;	//No optimization if we don't have 2 fanin signals each
	}
	
	
	signal a,b,c,d;	//a in sigVec1, c in sigVec2, b in both always (d the occurence of b in sigVec2)
	
	if(ntk.node_to_index(ntk.get_node(sigVec1.at(0))) == ntk.node_to_index(ntk.get_node(sigVec2.at(0)))){	//node 1 == node 3
		a = sigVec1.at(1);
		b = sigVec1.at(0);
		c = sigVec2.at(1);
		d = sigVec2.at(0);
	}else if(ntk.node_to_index(ntk.get_node(sigVec1.at(0))) == ntk.node_to_index(ntk.get_node(sigVec2.at(1)))){	//node 1 == node 4
		a = sigVec1.at(1);
		b = sigVec1.at(0);
		c = sigVec2.at(0);
		d = sigVec2.at(1);
	}else if(ntk.node_to_index(ntk.get_node(sigVec1.at(1))) == ntk.node_to_index(ntk.get_node(sigVec2.at(0)))){	//node 2 == node 3
		a = sigVec1.at(0);
		b = sigVec1.at(1);
		c = sigVec2.at(1);
		d = sigVec2.at(0);
	}else if(ntk.node_to_index(ntk.get_node(sigVec1.at(1))) == ntk.node_to_index(ntk.get_node(sigVec2.at(1)))){	//node 2 == node 4
		a = sigVec1.at(0);
		b = sigVec1.at(1);
		c = sigVec2.at(0);
		d = sigVec2.at(1);
	}else{
		return false;	//No distributivity possible
	}

	if(ntk.is_complemented(b) != ntk.is_complemented(d)){
		return false;	//Distributivity requires same complement in both occurences of b
	}

	if(ntk.is_complemented(sigVec.at(0)) && ntk.is_complemented(sigVec.at(1))){	//Annoying case with output inversion
		a = !a;
		c = !c;
	
		signal newSigVec2 = ntk.create_nand(a, c);	//Creates (ac) gate with correct complement

		signal nSig = ntk.create_and(b, newSigVec2);
		
		ntk.substitute_node(n, !nSig);	//replace current node with new signal

	}else{	//Easier case with only inversions on primary inputs
		if(ntk.is_complemented(sigVec.at(0))){	//Complement propagation
			a = !a;
		}
		if(ntk.is_complemented(sigVec.at(1))){
			c = !c;
		}
		
		signal newSigVec2 = ntk.create_and(a, c);	//Creates (ac) gate with correct complement
		
		ntk.substitute_node(ntk.get_node(sigVec.at(0)), b);	//b replaces fanin 1 of n
		ntk.substitute_node(ntk.get_node(sigVec.at(1)), newSigVec2);	//ac replaces fanin 2
	}

 	
	return true;
  }
  
  /* Try the three layer distributivity rule on node n. Return true if the network is updated. */
  bool try_3_layer_distributivity( node n )
  {
  	std::vector<signal> sigVec;
  	
	ntk.foreach_fanin( n, [&]( signal sig ){ sigVec.push_back(sig); });	//Extracts the fanin signals
	
	if(sigVec.size() != 2){
		return false;	//No optimization if we don't have 2 fanin signals or correct shape
	}

	if(ntk.level(ntk.get_node(sigVec.at(0)))+1 < ntk.level(ntk.get_node(sigVec.at(1)))){
		std::swap(sigVec.at(0), sigVec.at(1));
	}else if(ntk.level(ntk.get_node(sigVec.at(1)))+1 < ntk.level(ntk.get_node(sigVec.at(0)))){
		return false; //Nothing to gain from optimization
	}


	if(!ntk.is_complemented(sigVec.at(0))){
		return false;	//Wrong configuration for simplification
	}

  	std::vector<signal> sigVec2;
  	
	ntk.foreach_fanin( ntk.get_node(sigVec.at(0)), [&]( signal sig ){ sigVec2.push_back(sig); });	//Extracts the fanin signals
	
	if(sigVec2.size() != 2){
		return false;	//No optimization if we don't have 2 fanin signals or correct shape
	}
	
	if(ntk.level(ntk.get_node(sigVec2.at(0))) < ntk.level(ntk.get_node(sigVec2.at(1)))){
		std::swap(sigVec2.at(0), sigVec2.at(1));
	}

	if(ntk.is_complemented(sigVec.at(0))){
		return false;	//Wrong configuration for simplification
	}

  	std::vector<signal> sigVec3;
  	
	ntk.foreach_fanin( ntk.get_node(sigVec2.at(0)), [&]( signal sig ){ sigVec3.push_back(sig); });	//Extracts the fanin signals
	
	if(sigVec3.size() != 2){
		return false;	//No optimization if we don't have 2 fanin signals or correct shape
	}
	
	if(ntk.level(ntk.get_node(sigVec3.at(0))) < ntk.level(ntk.get_node(sigVec3.at(1)))){
		std::swap(sigVec3.at(0), sigVec3.at(1));
	}
	
	signal x2x4 = ntk.create_and(sigVec3.at(1), sigVec.at(1));
	signal gx2x4 = ntk.create_nand(sigVec3.at(0), x2x4);
	signal x3x4 = ntk.create_nand(!(sigVec2.at(1)), sigVec.at(1));
	signal newOutSig = ntk.create_nand(gx2x4, x3x4);
	
	ntk.substitute_node(n, newOutSig);
	
	return true;
  }

private:
  Ntk& ntk;
};

} // namespace detail

/* Entry point for users to call */
template<class Ntk>
void aig_algebraic_rewriting( Ntk& ntk )
{
  static_assert( std::is_same_v<typename Ntk::base_type, aig_network>, "Ntk is not an AIG" );

  depth_view dntk{ntk};
  detail::aig_algebraic_rewriting_impl p( dntk );
  p.run();
}

} /* namespace mockturtle */
