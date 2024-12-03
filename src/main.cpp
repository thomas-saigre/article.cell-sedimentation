/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

const double G_CST = 9.80665;

std::string computeGravityExpression(double time, double T_wait=0.3, double T_rotation = 3 /*, position initial_pos = SUPINE, position final_pos = PRONE*/)
{
    std::string gravity_expr;
    if (time >= T_rotation + T_wait)
        // gravity_expr = "{0,0,-" + G_CST + "}";
        gravity_expr = fmt::format("{{-{},0,0}}", G_CST);
    else if (time <= T_wait)
        gravity_expr = fmt::format("{{{},0,0}}", G_CST);
    else
    {
        double t_normalized = (time - T_wait) / (T_rotation);
        double angle_x = cos(M_PI * t_normalized),
               angle_z = sin(M_PI * t_normalized);
        // gravity_expr = "{" + std::to_string(G_CST*angle_x) + ",0," + std::to_string(G_CST*angle_z) + "}";
        gravity_expr = fmt::format("{{{},0,{}}}", G_CST * angle_x, G_CST * angle_z);
    }
    return gravity_expr;
}

template <int nDim, int OrderT, int OrderV, int OrderP>
int
runApplicationHeatFluidTwoToolboxes()
{
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<nDim,1>,
                                Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_heat_type;
    typedef FeelModels::FluidMechanics< Simplex<nDim,1>,
                                        Lagrange<OrderV, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderP, Scalar,Continuous,PointSetFekete> > model_fluid_type;
    typedef FeelModels::HeatFluid< model_heat_type,model_fluid_type> model_type;
    std::shared_ptr<model_type> heatFluidSupine( new model_type("heat-fluid") ),
                                heatFluid( new model_type("heat-fluid") );

    heatFluidSupine->init();
    heatFluidSupine->setStationary( true );
    heatFluidSupine->fluidModel()->setStationary( true );
    heatFluidSupine->heatModel()->setStationary( true );
    heatFluidSupine->printAndSaveInfo();

    std::ofstream wss_file("wss.csv");

    if (Environment::isMasterRank())
    {
        wss_file << "time,domain,AqueousHumor_Cornea,AqueousHumor_Iris\n";

        std::cout << "============================================================\n";
        std::cout << "Initial condition : Supine position (t = " << heatFluid->time() << ")\n";
        std::cout << "============================================================\n";
    }
    heatFluidSupine->solve();
    heatFluidSupine->exportResults();

    if (Environment::isMasterRank())
        wss_file << heatFluid->time();

    for ( std::string const& name : { "", "AqueousHumor_Cornea", "AqueousHumor_Iris"} )
    {
        auto wss = heatFluidSupine->fluidModel()->computeWallShearStress( name );
        auto wss_mean = mean( _range = elements(heatFluidSupine->fluidModel()->functionSpaceVelocity()->template meshSupport<0>()),
                              _expr = sqrt(inner(idv(wss), idv(wss))) );

        if (Environment::isMasterRank())
            wss_file << "," << wss_mean;
    }
    wss_file << std::endl;

    // heatFluid->updateTimeStep();

    std::string gravity_expr;

    auto const& u = heatFluidSupine->fluidModel()->fieldVelocity();
    auto const& p = heatFluidSupine->fluidModel()->fieldPressure();
    auto const& T = heatFluidSupine->heatModel()->fieldTemperature();

    heatFluid->init();
    *heatFluid->fluidModel()->fieldVelocityPtr() = u;
    *heatFluid->fluidModel()->fieldPressurePtr() = p;
    *heatFluid->heatModel()->fieldTemperaturePtr() = T;


    for ( heatFluid->startTimeStep() ; !heatFluid->timeStepBase()->isFinished(); heatFluid->updateTimeStep() )
    {
        if (heatFluid->worldComm().isMasterRank())
        {
            std::cout << "============================================================\n";
            std::cout << "time simulation: " << heatFluid->time() << "s\n";
            std::cout << "============================================================\n";
        }

        gravity_expr = computeGravityExpression(heatFluid->time());
        std::cout << "Gravity expression: " << gravity_expr << std::endl;
        heatFluid->updateGravityForce(gravity_expr);

        heatFluid->solve();
        heatFluid->exportResults();

        if (Environment::isMasterRank())
            wss_file << heatFluid->time();

        for ( std::string const& name : { "", "AqueousHumor_Cornea", "AqueousHumor_Iris"} )
        {
            auto wss = heatFluid->fluidModel()->computeWallShearStress( name );
            auto wss_mean = mean( _range = elements(heatFluid->fluidModel()->functionSpaceVelocity()->template meshSupport<0>()),
                                _expr = sqrt(inner(idv(wss), idv(wss))) );

            if (Environment::isMasterRank())
            {
                wss_file << "," << wss_mean;
            }
        }
    wss_file << std::endl;
    }
    wss_file.close();

    return !heatFluid->checkResults();
}


template <int nDim, int OrderT, int OrderV, int OrderP>
int
runApplicationHeatFluid()
{
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<nDim,1>,
                                Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_heat_type;
    typedef FeelModels::FluidMechanics< Simplex<nDim,1>,
                                        Lagrange<OrderV, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderP, Scalar,Continuous,PointSetFekete> > model_fluid_type;
    typedef FeelModels::HeatFluid< model_heat_type,model_fluid_type> model_type;
    std::shared_ptr<model_type> heatFluid( new model_type("heat-fluid") );
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;

    heatFluid->init();
    heatFluid->printAndSaveInfo();

    if ( !heatFluid->doRestart() )
        heatFluid->exportResults();


    std::ofstream wss_file("wss.csv");

    if (Environment::isMasterRank())
    {
        wss_file << "time,domain,AqueousHumor_Cornea,AqueousHumor_Iris\n";
    }

    for ( std::string const& name : { "", "AqueousHumor_Cornea", "AqueousHumor_Iris"} )
    {
        auto wss = heatFluid->fluidModel()->computeWallShearStress( name );
        auto wss_mean = mean( _range = elements(heatFluid->fluidModel()->functionSpaceVelocity()->template meshSupport<0>()),
                            _expr = sqrt(inner(idv(wss), idv(wss))) );

        if (Environment::isMasterRank())
        {
            wss_file << "," << wss_mean;
        }
    }

    // heatFluid->updateTimeStep();

    std::string gravity_expr;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;


    for ( heatFluid->startTimeStep() ; !heatFluid->timeStepBase()->isFinished(); heatFluid->updateTimeStep() )
    {
        if (heatFluid->worldComm().isMasterRank())
        {
            std::cout << "============================================================\n";
            std::cout << "time simulation: " << heatFluid->time() << "s\n";
            std::cout << "============================================================\n";
        }

        gravity_expr = computeGravityExpression(heatFluid->time());
        std::cout << "Gravity expression: " << gravity_expr << std::endl;
        heatFluid->updateGravityForce(gravity_expr);

        heatFluid->solve();
        heatFluid->exportResults();

        if (Environment::isMasterRank())
            wss_file << heatFluid->time();

        for ( std::string const& name : { "", "AqueousHumor_Cornea", "AqueousHumor_Iris"} )
        {
            auto wss = heatFluid->fluidModel()->computeWallShearStress( name );
            auto wss_mean = mean( _range = elements(heatFluid->fluidModel()->functionSpaceVelocity()->template meshSupport<0>()),
                                _expr = sqrt(inner(idv(wss), idv(wss))) );

            if (Environment::isMasterRank())
            {
                wss_file << "," << wss_mean;
            }
        }
    wss_file << std::endl;
    }
    wss_file.close();

    return !heatFluid->checkResults();
}

int
main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description heatfluidoptions( "application heat-fluid options" );
    heatfluidoptions.add( toolboxes_options("heat-fluid") );
    heatfluidoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1-P2P1" ), "discretization : P1-P2P1")
     ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=heatfluidoptions,
                     _about=about(_name="feelpp_toolbox_heatfluid_wss",
                                  _author="Thomas Saigre",
                                  _email="thomas.saigre@cemosis.fr"));

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    std::cout << __FILE__ << ":" << __LINE__ << std::endl;

    auto dimt = hana::make_tuple(hana::int_c<3>);

    auto discretizationt = hana::make_tuple( hana::make_tuple("P1-P2P1", hana::make_tuple( hana::int_c<1>,hana::int_c<2>,hana::int_c<1>) ));

    int status = 0;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension,&status]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _torder = std::decay_t<decltype(hana::at_c<0>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        constexpr int _uorder = std::decay_t<decltype(hana::at_c<1>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        constexpr int _porder = std::decay_t<decltype(hana::at_c<2>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        if ( dimension == _dim && discretization == _discretization )
                            status = runApplicationHeatFluid<_dim,_torder,_uorder,_porder>();
                    } );

    return status;
}
