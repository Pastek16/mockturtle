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
	
		if(ntk.is_on_critical_path(n) == 0){
			return false;	//Skip
		}
		
		std::vector<signal> sigVec;
		
		ntk.foreach_fanin( n, [&]( signal sig ){ sigVec.push_back(sig); });	//Extracts the fanin signals
		
		if(sigVec.size() != 2){
			return false;	//No optimization if we don't have 2 fanin signals
		}
		
		if(ntk.is_pi(ntk.get_node(sigVec.at(0))) && ntk.is_pi(ntk.get_node(sigVec.at(1)))){
			return false;
		}
		
		if(ntk.level(ntk.get_node(sigVec.at(0))) >= (ntk.level(ntk.get_node(sigVec.at(1)))+2) && ntk.is_complemented(sigVec.at(0)) == 0){
			std::swap(sigVec.at(0), sigVec.at(1));	//signal 1 being split
		}else if(ntk.level(ntk.get_node(sigVec.at(1))) >= (ntk.level(ntk.get_node(sigVec.at(0)))+2) && ntk.is_complemented(sigVec.at(1)) == 0){

		}else{
			return false;
		}
		

		//	Fanins of fanin
		std::vector<signal> sigVecBot;
		
		ntk.foreach_fanin( ntk.get_node(sigVec.at(1)), [&]( signal sig ){ sigVecBot.push_back(sig); });	//Extracts the fanin signals

		if(sigVecBot.size() != 2){
			return false;	//No optimization if we don't have 2 fanin signals
		}

		if(ntk.is_on_critical_path(ntk.get_node(sigVecBot.at(0))) != 0 && ntk.is_on_critical_path(ntk.get_node(sigVecBot.at(1))) == 0){
			std::swap(sigVecBot.at(0), sigVecBot.at(1));	//signal 1 being exchanged
		}else if(ntk.is_on_critical_path(ntk.get_node(sigVecBot.at(1))) != 0 && ntk.is_on_critical_path(ntk.get_node(sigVecBot.at(0))) == 0){

		}else{
			return false;
		}
		
		signal newOut = ntk.create_and(ntk.create_and(sigVec.at(0), sigVecBot.at(0)), sigVecBot.at(1));
		ntk.substitute_node(n, newOut);

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
	
	if(ntk.node_to_index(ntk.get_node(sigVec1.at(0))) == ntk.node_to_index(ntk.get_node(sigVec2.at(0)))){	//node 1 == node 3

	}else if(ntk.node_to_index(ntk.get_node(sigVec1.at(0))) == ntk.node_to_index(ntk.get_node(sigVec2.at(1)))){	//node 1 == node 4
		std::swap(sigVec2.at(0), sigVec2.at(1));
	}else if(ntk.node_to_index(ntk.get_node(sigVec1.at(1))) == ntk.node_to_index(ntk.get_node(sigVec2.at(0)))){	//node 2 == node 3
		std::swap(sigVec1.at(0), sigVec1.at(1));
	}else if(ntk.node_to_index(ntk.get_node(sigVec1.at(1))) == ntk.node_to_index(ntk.get_node(sigVec2.at(1)))){	//node 2 == node 4
		std::swap(sigVec1.at(0), sigVec1.at(1));
		std::swap(sigVec2.at(0), sigVec2.at(1));	//All swaps to ensure node 1 == node 3
	}else{
		return false;	//No distributivity possible
	}

	if(ntk.is_complemented(sigVec1.at(0)) != ntk.is_complemented(sigVec2.at(0))){
		return false;	//Distributivity requires same complement in both occurences of b
	}

	if(ntk.is_complemented(sigVec.at(0)) && ntk.is_complemented(sigVec.at(1))){	//Annoying case with output inversion
		sigVec1.at(1) = !sigVec1.at(1);
		sigVec2.at(1) = !sigVec2.at(1);
	
		signal newSigVec2 = ntk.create_nand(sigVec1.at(1), sigVec2.at(1));	//Creates (ac) gate with correct complement

		signal nSig = ntk.create_and(sigVec1.at(0), newSigVec2);
		
		ntk.substitute_node(n, !nSig);	//replace current node with new signal

	}else{	//Easier case with only inversions on primary inputs
		if(ntk.is_complemented(sigVec.at(0))){	//Complement propagation
			sigVec1.at(1) = !sigVec1.at(1);
		}
		if(ntk.is_complemented(sigVec.at(1))){
			sigVec2.at(1) = !sigVec2.at(1);
		}
		
		signal newSigVec2 = ntk.create_and(sigVec1.at(1), sigVec2.at(1));	//Creates (ac) gate with correct complement
		
		signal nSig = ntk.create_and(sigVec1.at(0), newSigVec2);
		ntk.substitute_node(n, nSig);	//ac replaces fanin 2
	}

 	
	return true;
  }
  
  /* Try the three layer distributivity rule on node n. Return true if the network is updated. */
  bool try_3_layer_distributivity( node n )
  {

	std::vector <signal> sigTop;
	
	ntk.foreach_fanin(n, [&](signal in){sigTop.push_back(in);});
	
	if(ntk.level(ntk.get_node(sigTop.at(0))) < ntk.level(ntk.get_node(sigTop.at(1)))){	//Deepest always on left
		std::swap(sigTop.at(0), sigTop.at(1));
	}
	if(ntk.level(ntk.get_node(sigTop.at(1)))+3 > ntk.level(ntk.get_node(sigTop.at(0)))){
		return false;
	}

	if(ntk.is_on_critical_path(ntk.get_node(sigTop.at(0))) == 0 || ntk.is_on_critical_path(ntk.get_node(sigTop.at(1))) != 0){
		return false;	//Wrong conditions
	}

	if(ntk.is_complemented(sigTop.at(0)) == 0){
		return false;	//Wrong conditions
	}
	
	std::vector <signal> sigMid;
	ntk.foreach_fanin(ntk.get_node(sigTop.at(0)), [&](signal in) {	sigMid.push_back(in);});

	if(ntk.level(ntk.get_node(sigMid.at(0))) < ntk.level(ntk.get_node(sigMid.at(1)))){	//Deepest always on left
		std::swap(sigMid.at(0), sigMid.at(1));
	}

	if(ntk.is_on_critical_path(ntk.get_node(sigMid.at(0))) == 0 || ntk.is_on_critical_path(ntk.get_node(sigMid.at(1))) != 0){
		return false;	//Wrong conditions
	}

	if(ntk.is_complemented(sigMid.at(0)) == 0){
		return false;	//Wrong conditions
	}
	
	std::vector <signal> sigBot;
	ntk.foreach_fanin(ntk.get_node(sigMid.at(0)), [&](signal in) {sigBot.push_back(in);});

	if(ntk.level(ntk.get_node(sigBot.at(0))) < ntk.level(ntk.get_node(sigBot.at(1)))){	//Deepest always on left
		std::swap(sigBot.at(0), sigBot.at(1));
	}

	if(ntk.is_on_critical_path(ntk.get_node(sigBot.at(0))) == 0 || ntk.is_on_critical_path(ntk.get_node(sigBot.at(1))) != 0){
		return false;
	}

	signal newOut = ntk.create_nand(ntk.create_nand(sigBot.at(0), ntk.create_and(sigBot.at(1), sigTop.at(1))) , ntk.create_nand(!sigMid.at(1), sigTop.at(1)));
	
	ntk.substitute_node(n, newOut);

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
