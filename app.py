from flask import Flask, render_template, request
import ref_water
import sediment_type
import numpy as np
import matplotlib

matplotlib.use('Agg')  # Use Agg backend
import matplotlib.pyplot as plt
from io import BytesIO
import base64

app = Flask(__name__)

def generate_plot(diameters, settling_velocities):
    plt.figure(figsize=(8, 6))
    plt.plot(diameters, settling_velocities, marker='o', linestyle='-')
    plt.xlabel('Particle Diameter (mm)')
    plt.ylabel('Settling Velocity (m/s)')
    plt.title('Settling Velocity vs Particle Diameter')
    plt.grid(True)

    # Convert plot to base64 encoded image
    img = BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()

    plt.close()

    return plot_url

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    # Get input from the form
    min_diameter = float(request.form['min_diameter']) / 1000.0
    max_diameter = float(request.form['max_diameter']) / 1000.0
    diameter_interval = float(request.form['diameter_interval']) / 1000.0
    particle_density = float(request.form['particle_density'])
    temperature = float(request.form['temperature'])
    salinity = float(request.form['salinity'])
    pH = float(request.form['pH'])

    water = ref_water.water(temperature=temperature, salinity=salinity, pH=pH)

    # Calculate sediment settling velocities for the given diameter range
    diameters_in_meter = [min_diameter + i * diameter_interval for i in range(int((max_diameter - min_diameter) / diameter_interval) + 1)]

    settling_velocities = [sediment_type.sand(water, d, particle_density).ws[0] for d in diameters_in_meter]
    settling_velocities = [abs(velocity) for velocity in settling_velocities]

    diameters = [d * 1000.0 for d in diameters_in_meter]  # Convert diameters back to millimeters

    # Generate plot
    plot_url = generate_plot(diameters, settling_velocities)

    return render_template('index.html',
                           diameters=diameters,
                           settling_velocities=settling_velocities,
                           plot_url=plot_url,
                           min_diameter=min_diameter*1000,
                           max_diameter=max_diameter*1000,
                           diameter_interval=diameter_interval*1000,
                           particle_density=particle_density,
                           temperature=temperature,
                           salinity=salinity,
                           pH=pH)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 8080))
    app.run(debug=False, host='0.0.0.0', port=port)
