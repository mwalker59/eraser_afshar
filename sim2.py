import sys
import json
import argparse
import numpy as np

from optical_system2 import OpticalSystem
from ray_graph import RayGraph


def main():
    # arguments

    parser = argparse.ArgumentParser(description="make a graph")

    parser.add_argument("sysname", help="name of optical system from JSON")
    parser.add_argument("--save", "-s", metavar='FILENAME', type=str)

    args = parser.parse_args()


    # load config file

    try:
        with open("optical_configs.json", "r") as f:
            configs = json.load(f)

        if args.sysname in configs:
            optical = configs[args.sysname]
        else:
            print("Config not found")
            sys.exit(0)

    except FileNotFoundError:
        print("Config file not found")
        sys.exit(1)


    # create image space

    rg = RayGraph(optical["graph_width"], optical["graph_height"],
                  0, optical["graph_height"]/2,
                  optical["field_width"], optical["field_height"])


    # create optical system object

    osys = OpticalSystem(qwp_upper=optical["qwp_upper"], qwp_lower=optical["qwp_lower"], wavelength=optical["wavelength"],
                         x_lens=optical["x_lens"], r_lens=optical["r_lens"], n_lens=optical["n_lens"],
                         x_pos = 2.0, y_delta=optical["y_delta_source"],
                         field_width=optical["field_width"], field_height=optical["field_height"],
                         scale_factor=optical["angle_scale"])


    # create sensor/block objects

    sensor_array = []

    for sensor in optical["sensors"]:
        # create object

        if "pol_axis" in sensor:
            sensor_array.append(Sensor(sensor["x"], sensor["y"], sensor["width"], sensor["name"], np.array(sensor["pol_axis"])))
        else:
            sensor_array.append(Sensor(sensor["x"], sensor["y"], sensor["width"], sensor["name"]))


    # render polarizations

    def rotate_pol(base_vector, theta_rad):
        rotation_matrix = np.array([
            [np.cos(theta_rad), -np.sin(theta_rad)],
            [np.sin(theta_rad),  np.cos(theta_rad)]
        ])

        return rotation_matrix @ base_vector


    for pol in optical["polarizations"]:
        if pol["pol"] == "H": pol2 = OpticalSystem.pol_H
        if pol["pol"] == "V": pol2 = OpticalSystem.pol_V
        if pol["pol"] == "D": pol2 = OpticalSystem.pol_D
        if pol["pol"] == "A": pol2 = OpticalSystem.pol_A

        color = [pol["color_vec"][0], pol["color_vec"][1], pol["color_vec"][2]]

        pol2 = rotate_pol(pol2, optical["polarization_rotation"] * (2*np.pi / 360.0))
        render_interference_image(osys, pol2, optical, rg, color, sensor_array)


    # draw sensors

    for sensor in optical["sensors"]:
        # print sensor on image

        rg.draw_line(sensor["x"], sensor["y"] + sensor["width"], sensor["x"], sensor["y"] - sensor["width"], [1.0, 1.0, 1.0])


    for sensor in sensor_array:
        print(sensor.get_name() + ": " + str(sensor.get_accum_power()))


    if args.save is not None:
        rg.save(sensor_array[0].get_accum_power(), sensor_array[1].get_accum_power(), args.save)

    else:
        rg.show()


class Sensor:
    def __init__(self, x, y, width, name, polarizer_axis=None):
        self.x = x
        self.y = y
        self.width = width
        self.name = name
        self.pol_axis = polarizer_axis
        self.accum_power = 0.0


    def intersect(self, points, point_idx, inten, polarization):
        if points[point_idx][0] > self.x or points[point_idx + 1][0] < self.x:
            return (False, 0, 0)

        x_delta = points[point_idx + 1][0] - points[point_idx][0]
        y_delta = points[point_idx + 1][1] - points[point_idx][1]

        x_frac = (self.x - points[point_idx][0]) / x_delta
        y_intersect = points[point_idx][1] + y_delta * x_frac

        if y_intersect < self.y + self.width and y_intersect > self.y - self.width:
            points[point_idx + 1][0] = self.x
            points[point_idx + 1][1] = y_intersect

            if self.pol_axis is not None:
                dot = np.dot(polarization, self.pol_axis)
                self.accum_power += inten * (dot**2).real
                print("Add to " + self.name + ": " + str(inten * (dot**2).real), self.accum_power)
            else:
                self.accum_power += inten

            return (True, self.x, y_intersect)

        else:
            return (False, 0, 0)


    def get_name(self):
        return self.name

    def get_accum_power(self):
        return self.accum_power


def render_interference_image(osys, polarization, optical, rg, color_vector, sensor_array):
    # private function to walk through angles and draw

    def walk_angles(high_low):
        angle = optical["low_angle"]

        while angle <= optical["high_angle"]:
            # trace ray

            points, inten = osys.trace_source(polarization, high_low, angle)


            # draw segments

            for point_idx in range(0, len(points) - 1):
                # check whether line intersects a sensor

                abort = False

                for sensor in sensor_array:
                    abort_temp, new_x, new_y = sensor.intersect(points, point_idx, inten, polarization)
                    if abort_temp: abort = True


                # get start and end, draw line

                x0 = points[point_idx][0]
                y0 = points[point_idx][1]

                x1 = points[point_idx + 1][0]
                y1 = points[point_idx + 1][1]

                rg.draw_line(x0, y0, x1, y1, inten * np.array(color_vector))


                # check whether segment intersected a sensor/block

                if abort:
                    break


            # advance to next angle

            angle += optical["angle_step"]


    # draw rays from sources

    if optical["upper"] > 0:
        walk_angles("upper")

    if optical["lower"] > 0:
        walk_angles("lower")


if __name__ == "__main__":
    main()
