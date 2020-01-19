/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "juiceIntegratorAnalysisTools.h"

std::string getCurrentRootPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string currentPath = filePath_.substr( 0, filePath_.length( ) -
                                                ( std::string( "juiceIntegratorAnalysisTools.cpp" ).length( ) ) );

    return currentPath;
}

double getClosestApproachTime( std::map< double, Eigen::VectorXd > numericalSolution )
{
    //! STUDENT CODE TASK: compute epoch of closest approach


    double closestApproachTime =0.0;
    double mindistance = std::numeric_limits< double >::infinity( );

    for(auto Iterator : numericalSolution)
        {
            double currentTime = Iterator.first;
            Eigen::Vector6d currentstate=Iterator.second;
            //finding absolute distance
            double currentdistance = currentstate.segment(0,3).norm();

            if (currentdistance < mindistance)
                {
                closestApproachTime = currentTime;
            }
        }
    return closestApproachTime ;
}

std::shared_ptr< IntegratorSettings< double > > getFixedStepSizeIntegratorSettings(
        const double initialTime, const double timeStep )
{
     //! STUDENT CODE TASK: create fixed step RKF7(8) integrator settings (uncomment next line, and fill in .... ).
    return std::make_shared< RungeKuttaVariableStepSizeSettings< double>> (initialTime, timeStep, RungeKuttaCoefficients::rungeKuttaFehlberg78, timeStep, timeStep, 1,1);
}

AccelerationMap getUnperturbedAccelerations(
        const std::string& centralBody, const NamedBodyMap& bodyMap )
{
    // Define list of acceleration settings acting on JUICE
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfJuice;
     //! STUDENT CODE TASK: define acceleration settings for unperturbed dynamics (fill up accelerationsOfJuice)
    // Create and return full list of accelerations
    accelerationsOfJuice[ centralBody ].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::central_gravity ) );
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "JUICE" ] = accelerationsOfJuice;
    return createAccelerationModelsMap(
                bodyMap, accelerationSettings, { "JUICE" }, { centralBody } );
}

AccelerationMap getPerturbedAccelerations(
        const std::string& centralBody, const NamedBodyMap& bodyMap )
{
    // Define list of acceleration settings acting on JUICE
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfJuice;

    //! STUDENT CODE TASK: define acceleration settings for perturbed dynamics (fill up accelerationsOfJuice)
    accelerationsOfJuice[ centralBody ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfJuice[ "Jupiter" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 6, 0 ) );
    accelerationsOfJuice[ "Io" ].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::central_gravity ) );
    accelerationsOfJuice[ "Europa" ].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::central_gravity ) );
    accelerationsOfJuice[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::central_gravity ) );
    accelerationsOfJuice[ "Saturn" ].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::central_gravity ) );

    accelerationsOfJuice[ "Sun"].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::cannon_ball_radiation_pressure ) );
    if(centralBody == "Ganymede" )
    {
    accelerationsOfJuice[ "Ganymede"].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::aerodynamic ) );
    accelerationsOfJuice[ "Callisto"].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::central_gravity ) );
    }
    if(centralBody == "Callisto" )
    {
    accelerationsOfJuice[ "Ganymede"].push_back( std::make_shared< AccelerationSettings >(basic_astrodynamics::central_gravity ) );
    }

    // Create and return full list of accelerations
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "JUICE" ] = accelerationsOfJuice;
    return createAccelerationModelsMap(
                bodyMap, accelerationSettings, { "JUICE" }, { centralBody } );
}

void writePropagationResultsToFile(
        std::shared_ptr< SingleArcDynamicsSimulator< > >& dynamicsSimulator,
        const std::string& fileOutputIdentifier )
{
    std::string outputFolder = getCurrentRootPath( ) + "SimulationOutput";

    input_output::writeDataMapToTextFile(
                dynamicsSimulator->getEquationsOfMotionNumericalSolution( ), fileOutputIdentifier + "_numerical_states.dat", outputFolder );
    input_output::writeDataMapToTextFile(
                dynamicsSimulator->getDependentVariableHistory( ), fileOutputIdentifier + "_dependent.dat", outputFolder );

}


void writePropagationResultsAndBenchmarkDifferenceToFile(
        std::shared_ptr< SingleArcDynamicsSimulator< > >& dynamicsSimulator,
        const std::string& fileOutputIdentifier,
        const std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > benchmarkInterpolator )
{
    std::map< double, Eigen::VectorXd > numericalSolution = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableSolution =
            dynamicsSimulator->getDependentVariableHistory( );

    std::map< double, Eigen::VectorXd > benchmarkDifference;
    for( auto stateIterator : numericalSolution )
    {
        benchmarkDifference[ stateIterator.first ] =
                stateIterator.second - benchmarkInterpolator->interpolate( stateIterator.first );
    }

    std::string outputFolder = getCurrentRootPath( ) + "SimulationOutput";

    input_output::writeDataMapToTextFile(
                numericalSolution,
                fileOutputIdentifier + "_numerical_states.dat", outputFolder );
    input_output::writeDataMapToTextFile(
                benchmarkDifference, fileOutputIdentifier + "_benchmark_difference.dat", outputFolder );
    if( dependentVariableSolution.size( ) > 0 )
    {
        if( dependentVariableSolution.begin( )->second.rows( ) > 0 )
        {
            input_output::writeDataMapToTextFile(
                        dependentVariableSolution,
                        fileOutputIdentifier + "_dependent_states.dat", outputFolder );
        }
    }
}

std::map< double, Eigen::VectorXd > getDifferenceWrtKeplerOrbit(
        const std::map< double, Eigen::VectorXd >& numericalSolution,
        const double centralBodyGravitationalParameter )
{
    Eigen::Vector6d initialKeplerianElements =
            orbital_element_conversions::convertCartesianToKeplerianElements(
                Eigen::Vector6d( numericalSolution.begin( )->second ), centralBodyGravitationalParameter );
    double initialTime = numericalSolution.begin( )->first;

    std::map< double, Eigen::VectorXd > keplerianSolutionDifference;
    for( auto stateIterator : numericalSolution )
    {
        keplerianSolutionDifference[ stateIterator.first ] =
                orbital_element_conversions::convertKeplerianToCartesianElements(
                    propagateKeplerOrbit(
                        initialKeplerianElements, stateIterator.first - initialTime, centralBodyGravitationalParameter ),
                    centralBodyGravitationalParameter ) - stateIterator.second;
    }
    return keplerianSolutionDifference;
}


void writePropagationResultsAndAnalyticalSolutionToFile(
        std::shared_ptr< SingleArcDynamicsSimulator< > >& dynamicsSimulator,
        const std::string& fileOutputIdentifier,
        const double centralBodyGravitationalParameter )
{
    std::map< double, Eigen::VectorXd > numericalSolution = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableSolution =
            dynamicsSimulator->getDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > keplerianSolutionDifference =
            getDifferenceWrtKeplerOrbit( numericalSolution, centralBodyGravitationalParameter );

    std::string outputFolder = getCurrentRootPath( ) + "SimulationOutput";

    input_output::writeDataMapToTextFile(
                dynamicsSimulator->getEquationsOfMotionNumericalSolution( ),
                fileOutputIdentifier + "_numerical_states.dat", outputFolder );
    input_output::writeDataMapToTextFile(
                keplerianSolutionDifference, fileOutputIdentifier + "_keplerian_difference.dat", outputFolder );
    if( dependentVariableSolution.size( ) > 0 )
    {
        if( dependentVariableSolution.begin( )->second.rows( ) > 0 )
        {
            input_output::writeDataMapToTextFile(
                        dependentVariableSolution,
                        fileOutputIdentifier + "_dependent_states.dat", outputFolder );
        }
    }
}
