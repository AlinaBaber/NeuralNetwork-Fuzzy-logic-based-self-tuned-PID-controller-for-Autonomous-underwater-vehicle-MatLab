# Neural Network & Fuzzy Logic-based Self-tuned PID Controller for Autonomous Underwater Vehicle (AUV)

This project implements an advanced control system using a **Neural Network-Fuzzy Logic-based Self-tuned PID Controller** to optimize the performance and stability of an Autonomous Underwater Vehicle (AUV). This system combines neural networks, fuzzy logic, and a PID controller to adaptively control AUV movement in dynamic underwater environments.

## Project Overview

Underwater vehicles face challenging conditions due to unpredictable currents, variable pressure, and high resistance. Traditional PID controllers struggle with these non-linearities. This project proposes an intelligent, adaptive PID controller that uses:
- **Neural Networks** for predictive adjustments.
- **Fuzzy Logic** to tune the PID parameters in real time.
- **PID Control** to achieve precise positioning and stability for the AUV.

## System Components
![image](https://github.com/user-attachments/assets/eb1dcc42-40f4-4488-a086-8b616abbc1aa)

1. **Neural Network Model**:
   - Trains on environmental data to predict system behavior and adaptively adjust controller parameters.
   - Enables fast response to non-linear disturbances.

2. **Fuzzy Logic Controller**:
   - Uses a set of fuzzy rules to adjust PID parameters dynamically.
   - Provides real-time tuning by evaluating error and rate of change in error.

3. **PID Controller**:
   - Core feedback controller that maintains the desired trajectory by adjusting based on real-time error input.

## Project Features

- **Adaptive Control**: Real-time adjustments to PID parameters ensure stability and responsiveness in changing underwater conditions.
- **MATLAB Implementation**: Fully implemented in MATLAB with Simulink support for easy simulation and visualization.
- **Modular Design**: The controller's modular setup makes it flexible for various AUV models and operating conditions.
## Results:

## Prerequisites

- MATLAB (R2021a or later)
- Control System Toolbox
- Fuzzy Logic Toolbox
- Deep Learning Toolbox (for neural network training)

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/NeuralNetwork-Fuzzy-logic-based-self-tuned-PID-controller-for-Autonomous-underwater-vehicle-MatLab.git

License
This project is licensed under the MIT License - see the LICENSE file for details.
