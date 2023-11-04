#ifndef MINT2DHOOK_H
#define MINT2DHOOK_H

#include <mutex>
#include <thread>

#include "PhysicsHook.h"

#include "polyscope/polyscope.h"

#include "polyscope/surface_mesh.h"

#include <Eigen/Core>


#include "VizHelper.h"

#include "UtilsMisc.h"

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

// #include <fstream>
#include <sys/stat.h>

enum Field_View { vec_norms, delta_norms, vec_dirch, moment_dirch, primal_curl_residual, sym_curl_residual, gui_free, Element_COUNT };




class Mint2DHook : public virtual PhysicsHook
{
public:

    Mint2DHook() : PhysicsHook() {
    //   current_element = Field_View::vec_norms;
    }

    virtual void drawGUI()
    {
        //	ImGui::SliderFloat("k scale", &k_scale, 0.0001f, 2.0f, "k * %.3f");
        //	ImGui::SliderFloat("dt scale", &dt_scale, 0.0001f, 2.0f, "dt * %.3f");
            // ImGui::InputFloat("k scale", &k_scale);
            // ImGui::InputFloat("dt scale", &dt_scale);

        ImGui::InputDouble("Smoothness Weight", &w_smooth);
        ImGui::InputDouble("S Perp Weight", &w_s_perp);
        ImGui::InputDouble("Curl Weight", &w_curl);
        ImGui::InputDouble("Bound Weight", &w_bound);


/// Whatever maybe make this a dropdown eventually 
// From line 556 of imgui demo: https://skia.googlesource.com/external/github.com/ocornut/imgui/+/refs/tags/v1.73/imgui_demo.cpp
            const char* element_names[Field_View::Element_COUNT] = { "Vector Norms", "Delta Norms", "Vector Dirichlet", "Symmetric Dirichlet", "Vector Curl", "Symmetric Curl", "free" };
            const char* current_element_name = (current_element >= 0 && current_element < Field_View::Element_COUNT) ? element_names[current_element] : "Unknown";
            ImGui::PushItemWidth(300);
            ImGui::SliderInt("Shading Mode", (int *) &current_element, 0, Element_COUNT - 1, current_element_name);
            ImGui::PopItemWidth();
          
          
          
            // ImGui::SameLine(); // ImGui::HelpMarker("Using the format string parameter to display a name instead of the underlying integer.");

    }

    virtual void updateRenderGeometry()
    {



      vc.d().frames.resize(frames.rows(), 3);
      vc.d().frames << frames, Eigen::MatrixXd::Zero(frames.rows(), 1);
      vc.d().deltas = deltas;

      vc.updateVizState();

      vc.d().vec_curl = curls_primal;
      vc.d().sym_curl = curls_sym;

      vc.d().frame_smoothness = smoothness_primal;
      vc.d().moment_smoothness = smoothness_sym;



    }





// Need to fill out viewer for each of: Field_View { vec_dirch, moment_dirch, sym_curl_residual, primal_curl_residual,
    virtual void renderRenderGeometry()
    {
		polyscope::getSurfaceMesh()->updateVertexPositions(renderP);
        
        polyscope::getSurfaceMesh()->centerBoundingBox();
        polyscope::getSurfaceMesh()->resetTransform();

        switch (current_element) 
        { 
            case Field_View::vec_norms:
                { 
                  polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_norms", vc.d().frame_norms)->setEnabled(true);
                } 
                break;

            case Field_View::delta_norms:
                { 
                  polyscope::getSurfaceMesh()->addFaceScalarQuantity("delta_norms", vc.d().delta_norms)->setEnabled(true);
                } 
                break; 
              
            case Field_View::primal_curl_residual: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("primal_curl_residual", vc.d().vec_curl)->setEnabled(true);
                } 
                break;             
            case Field_View::sym_curl_residual: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("sym_curl_residual", vc.d().sym_curl)->setEnabled(true);
                } 
                break;           
            case Field_View::vec_dirch: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_dirch", vc.d().frame_smoothness)->setEnabled(true);
                } 
                break; 
            case Field_View::moment_dirch: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("moment_dirch", vc.d().moment_smoothness)->setEnabled(true);
                } 
                break; 

            default: 
                { 
                    // std::cout << "Unknown color!"; 
                } 
                break; 
        } 
      

  
        auto vectors = polyscope::getSurfaceMesh()->addFaceVectorQuantity("frames", vc.getFrames()); //   ( ((N.array()*0.5)+0.5).eval());
        vectors->vectorColor = glm::vec3(.7,.7,.7);
        
        // vectors->setVectorLengthScale(1., true);
        // vectors->setEnabled(true);
        // vectors->setVectorColor(glm::vec3(.7,.7,.7));
        
        polyscope::requestRedraw();   

        // std::ifstream file;
        std::string cur_log_file =  cur_log_folder + "/cur_file_iter_" + std::to_string(10000 + cur_iter) + ".png";

        struct stat buffer;   
        bool exists = (stat(cur_log_file.c_str(), &buffer) == 0); 

        if (!exists)
        polyscope::screenshot(cur_log_file, true);

        // // opening the file
        // file.open(cur_log_file.c_str(), 'r');

        // if (file == NULL)
              




    }




    virtual ~Mint2DHook()
    {      
    }

    public: 
            // static auto func;
        double w_bound;
        double w_smooth;
        double w_smooth_vector; 
        double w_curl;
        double w_s_perp;

        double w_attenuate;

        Field_View current_element;


    protected:
        VizHelper::VizCache vc;

        Eigen::MatrixXd V; // #V-by-3 3D vertex positions
        Eigen::MatrixXi F; // #F-by-3 indices into V
        Eigen::MatrixXd P; //  = tutte_embedding(V, F); // #V-by-2 2D vertex positions


        Eigen::MatrixXd renderP;
        Eigen::MatrixXi renderF;

        Eigen::MatrixXd frames;
        Eigen::MatrixXd deltas;

        Eigen::MatrixXd frames_orig;

        Eigen::VectorXd curls_sym;
        Eigen::VectorXd curls_primal; 
        Eigen::VectorXd smoothness_primal;
        Eigen::VectorXd smoothness_sym;

        std::string cur_mesh_name;
        std::string cur_log_folder;


        int max_iters = 5000;
        int cur_iter = 0;
        int inner_loop_iter = 0;
        double convergence_eps = 1e-10;
        double identity_weight = 1e-6;

        TinyAD::LinearSolver<double> solver;

          int buffer;

};
    
//     // virtual bool showSimButtons()
//     // {
//     //     return true;
//     // }

//     /*
//     * Runs when the user redraws/interacts with the GUI.
//     */
//     virtual void drawGUI() = 0;

//     /*
//      * Runs once when the simulation is initialized, and again each time the user resets the simulation.
//      */
//     virtual void initSimulation() = 0;

//     /*
//     * Called every once in a while (and always before every simulation step) even if the simulation is paused.
//     * Use this to update the visualization in response to user input etc.
//     */
//     virtual void tick() {};

//     /*
//      * Takes one simulation "step." You can do whatever you want here, but the granularity of your computation should 
//      * be small enough that the user can view/pause/kill the simulation at interactive rates.
//      * This method *must* be thread-safe with respect to renderRenderGeometry() (easiest is to not touch any rendering
//      * data structures at all).
//      */
//     virtual bool simulateOneStep() = 0;

//     /*
//      * Update the rendering data structures here. This method will be called in alternation with simulateOneStep().
//      * This method blocks rendering in the viewer, so do *not* do extensive computation here (leave it to 
//      * simulateOneStep()).
//      */
//     virtual void updateRenderGeometry() = 0;

//     /*
//      * Perform any actual rendering here. This method *must* be thread-safe with respect to simulateOneStep().
//      * This method runs in the same thread as the viewer and blocks user IO, so there really should not be any
//      * extensive computation here or the UI will lag/become unresponsive (the whole reason the simulation itself
//      * is in its own thread.)
//      */
//   //  virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) = 0;
//     virtual void renderRenderGeometry() = 0;

//     /*
//     * Called when the user clicks on the simulation panel.
//     * This method is called in the *rendering thread*. If you need to make changes to the simulation state, you
//     * should stash the mouse click in a message queue and deal with it in the simulation thread.
//     */
//  //   virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button) { return false; }

//     /*
//     * Called when the user unclicks the mouse on the simulation panel.
//     * This method is called in the *rendering thread*. If you need to make changes to the simulation state, you
//     * should stash the mouse click in a message queue and deal with it in the simulation thread.
//     */
//   //  virtual bool mouseReleased(igl::opengl::glfw::Viewer &viewer, int button) { return false; }

//     /*
//     * Called when the user drags the mouse on the simulation panel.
//     * This method is called in the *rendering thread*. If you need to make changes to the simulation state, you
//     * should stash the mouse click in a message queue and deal with it in the simulation thread.
//     */
//   //  virtual bool mouseMoved(igl::opengl::glfw::Viewer &viewer, int button) { return false; }

//     /*
//      * Runs the simulation, if it has been paused (or never started).
//      */
//     void run()
//     {
//         status_mutex.lock();
//         please_pause = false;
//         status_mutex.unlock();
//     }

//     /*
//      * Resets the simulation (and leaves it in a paused state; call run() to start it).
//      */
//     void reset()
//     {
//         killSimThread();
//         please_die = running = false;
//         please_pause = true;
//         initSimulation();

//         updateRenderGeometry();
//         sim_thread = new std::thread(&PhysicsHook::runSimThread, this);
//     }

//     /*
//      * Pause a running simulation. The simulation will pause at the end of its current "step"; this method will not
//      * interrupt simulateOneStep mid-processing.
//      */
//     void pause()
//     {
//         status_mutex.lock();
//         please_pause = true;
//         status_mutex.unlock();
//     }

//     bool isPaused()
//     {
//         bool ret = false;
//         status_mutex.lock();
//         if(running && please_pause)
//             ret = true;
//         status_mutex.unlock();
//         return ret;
//     }

//     void render()
//     {
//         render_mutex.lock();
//         renderRenderGeometry();
//         render_mutex.unlock();
//     }

//     bool showSimButtons()
//     {
//         return show_sim_buttons;
//     }
    
// protected:
//     void runSimThread()
//     {
//         status_mutex.lock();
//         running = true;
//         status_mutex.unlock();

//         bool done = false;
//         while (!done)
//         {
//             tick();

//             status_mutex.lock();
//             bool pausenow = please_pause;
//             status_mutex.unlock();
//             if (pausenow)
//             {
//                 std::this_thread::sleep_for(std::chrono::milliseconds(10));
//             }
//             else
//             {
//                 done = simulateOneStep();
//                 render_mutex.lock();
//                 updateRenderGeometry();
//                 render_mutex.unlock();
//             }
//             status_mutex.lock();
//             if (please_die)
//                 done = true;
//             status_mutex.unlock();
//         }

//         status_mutex.lock();
//         running = false;
//         status_mutex.unlock();
//     }

//     void killSimThread()
//     {
//         if (sim_thread)
//         {
//             status_mutex.lock();
//             please_die = true;
//             status_mutex.unlock();
//             sim_thread->join();
//             delete sim_thread;
//             sim_thread = NULL;
//         }
//     }

//     std::thread *sim_thread;
//     bool please_pause;
//     bool please_die;
//     bool running;
//     bool show_sim_buttons = true;
//     std::mutex render_mutex;
//     std::mutex status_mutex;
// };

#endif
