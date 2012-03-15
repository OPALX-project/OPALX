//
// C++ Implementation: ParticleLayoutFromGrid
//
// Description: 
//
//
//
// Author: Roman Geus <geus@maxwell>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
#include <cassert>
#include "ParticleLayoutFromGrid.h"

void ParticleLayoutFromGrid::update(ParticleBase< ParticleLayoutFromGrid >& particles) {
    unsigned num_procs = Ippl::getNodes();
    unsigned my_id = Ippl::myNode();
    size_t num_my_particles   = particles.getLocalNum();

    // Delete particles in destroy list, update local num
    size_t num_destroyed_particles = particles.getDestroyNum();
    particles.performDestroy();
    num_my_particles -= num_destroyed_particles;
    
    // FIXME: This should be done in ParticleBase::performDestroy()
    particles.setLocalNum(num_my_particles);

    // Apply boundary conditions to the particle positions
    // apply_bconds(particles.R);
    
    // Redistribute particles, such that particles that have moved outside the 
    // local domain are moved to the corresponding processor.
    unsigned num_active_procs = grid_->give_number_of_active_processes();
    if (num_active_procs > 1)
        num_my_particles = redistribute_particles(particles);
        
    // Update total number of particles in particle container
    size_t num_total_particles;
    reduce(num_my_particles, num_total_particles, OpAddAssign());
    particles.setTotalNum(num_total_particles);
}

void ParticleLayoutFromGrid::apply_bconds(ParticlePos_t& R)
{
    Inform msg("update");
    size_t num_my_particles = R.size();
    
    for (size_t ip = 0; ip < num_my_particles; ++ ip) {
        D3vector pos(R[ip](0), R[ip](1), R[ip](2));
        if (!geom_domain_->point_in_domain(pos)) {
            msg << "Particle " << ip << " moved outside domain." << endl;
        }
    }    
}

size_t ParticleLayoutFromGrid::redistribute_particles(ParticleBase< ParticleLayoutFromGrid >& particles) 
{
    unsigned num_procs = Ippl::getNodes();
    unsigned my_id = Ippl::myNode();
    size_t num_my_particles = particles.getLocalNum();
    
    // FIXME: How do I test if I am an active processor?
    // FIXME: Assume that active processors are assigned the ids 0..num_active_procs-1.
    unsigned num_active_procs = grid_->give_number_of_active_processes();

    // Non-active processors return here!
    if (my_id >= num_active_procs)
        return num_my_particles;
    
    std::vector<D3vector> bb_min;
    std::vector<D3vector> bb_max;
    for (int p = 0; p < num_active_procs; ++ p) {
        double xmin[3], xmax[3];
        grid_->localbox(xmin, xmax, p);
        bb_min.push_back(D3vector(xmin[0], xmin[1], xmin[2]));
        bb_max.push_back(D3vector(xmax[0], xmax[1], xmax[2]));
    }
    // FIXME: Bounding boxes could be sorted, such that neighbors are close the beginning.
    // FIXME: Should the local bounding box be removed from this list?
    
    // num_active_procs Message objects
    std::vector<Message> messages(num_active_procs);
    // List of particles for each message
    std::vector<size_t>* put_list = new std::vector<size_t>[num_active_procs];
    
    for (size_t ip = 0; ip < num_my_particles; ++ ip) {
        D3vector pos(particles.R[ip][0], particles.R[ip][1], particles.R[ip][2]);
        
        if (!is_local_pos(pos)) 
        {
            // Particle has moved outside the local domain
            
            // Identify processor which should receive the particle (linear search)
            unsigned target_id = num_active_procs;
            for (unsigned p = 0; p < num_active_procs; ++ p) {
                if (is_inside_box(bb_min[p], bb_max[p], pos)) {
                    target_id = p;
                    break;
                }
            }
            assert(target_id != num_active_procs && target_id != my_id);
            put_list[target_id].push_back(ip);
            
            // Mark particle for deletion
            particles.destroy(1, ip);
        }      
    }
    
    // Send particles to their destination nodes
    int tag = Ippl::Comm->next_tag(P_SPATIAL_TRANSFER_TAG, P_LAYOUT_CYCLE);
    for (int p = 0; p < num_active_procs; ++ p) {
        if (p != my_id) {
            // Put data for particles on this put list into message and add sentinel
            particles.putMessage(messages[p], put_list[p]);
            particles.putMessage(messages[p], (size_t) 0, (size_t) 0);
            // Send message
            // Note: The last parameter is set to false, to prevent send() from destroying the Message object
            Ippl::Comm->send(&messages[p], p, tag, false);
        }
        put_list[p].clear();
    }

    // Delete sent particles.
    num_my_particles -= particles.getDestroyNum();
    particles.performDestroy();
    
    // FIXME: This should be done in ParticleBase::performDestroy()
    particles.setLocalNum(num_my_particles);
    
    // Receive particles
    unsigned sendnum = num_active_procs - 1;
    while (sendnum-- > 0) {
        int node = Communicate::COMM_ANY_NODE;
        Message *recmsg = Ippl::Comm->receive_block(node, tag);
        size_t recvd;
        size_t total_recvd = 0;
        while ((recvd = particles.getMessage(*recmsg)) > 0)
            total_recvd += recvd;
        num_my_particles += total_recvd;
        delete recmsg;
    }
    
    // FIXME: This should be done in ParticleBase::getMessage()
    particles.setLocalNum(num_my_particles);
    
    delete[] put_list;
    return num_my_particles;
}
