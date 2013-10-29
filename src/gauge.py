#!/usr/bin/env python
"""

The Gauge widget draws a semi-circular gauge. You supply limits,
shaded regions, names and the current value, and invoke it like this:


    from pylab import figure, show
    
    current_value = -4.0
    limits = [-1.0,1.0,1,0.1]
    zone_colour = [[-1.0,0.0,'r'],[0.0,0.5,'y'],[0.5,1.0,'g']]
    attribute_name = "Rx MOS (24h)"
    
    
    graph_height = 1.6
    graph_width = 2.4
    fig_height = graph_height
    fig_width = graph_width
    
    fig = figure(figsize=(fig_width, fig_height ))
    
    rect = [(0.0/fig_width), (0.2/fig_height),
        (graph_width/fig_width), (graph_height/fig_height)]
    
    gauge = Gauge(fig, rect,
        xlim=( -0.1, graph_width+0.1 ),
        ylim=( -0.4, graph_height+0.1 ),
        xticks=[],
        yticks=[],
        )
    gauge.set_axis_off()
    fig.add_axes(gauge)
    
    show()

"""
from __future__ import division
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import math
import types

from math import pi


class Gauge(Axes):
    def __init__(self, *args, **kwargs):
        Axes.__init__(self, *args, **kwargs)

        #Perform Checking
        if( limits[0] == limits[1] ):
            raise ValueError('identical_limits_exception: %s'%limits)
        if( limits[1] > limits[0] ):
            graph_positive = True
        else: #Swap the limits around
            graph_positive = False
            limits[0], limits[1] = limits[1] = limits[0]
        
        #There must be an integer number of minor ticks for each major tick
        if not( ((limits[2]/limits[3]) % 1.0) * limits[3] == 0 ): 
            raise ValueError('bad_tick_spacing_exception')
        
        if( limits[2] <= 0 or
            limits[3] <= 0 or
            limits[2] < limits[3] or
            limits[3] > abs(limits[1]-limits[0]) ):
            raise ValueError('bad_limits_exception:%s' % limits)
        
        for zone in zone_colour:
            if( zone[0] > zone[1] ): #Swap the zones so zone[1] > zone[0]
                zone[0], zone[1] = zone[1], zone[0]
            if( zone[1] < limits[0] or zone[0] > limits[1] ):
                raise ValueError('bad_zone_exception'%zone)
            if( zone[0] < limits[0] ):
                zone[0] = limits[0]
            if( zone[1] > limits[1] ):
                zone[1] = limits[1]

        #Draw the arch
        for zone in zone_colour:
            self.draw_arch(limits, zone[0], zone[1], zone[2], False, 
                           graph_positive)
        self.draw_arch(limits, limits[0], limits[1], 'black', True, 
                       graph_positive)
        self.draw_ticks(limits, graph_positive)
        self.draw_needle(current_value, limits, graph_positive)
        self.draw_bounding_box(graph_width, graph_height)
        self.text(0.0, 0.2, attribute_name, size=10, va='bottom', 
                  ha='center')

        #The black dot.
        p = self.plot([0.0],[0.0],'.', color='#000000')


    def draw_arch(self, limits, start_angle, end_angle, colour, border, graph_positive):
        x_vect = []
        y_vect = []
        if( graph_positive ):
            start = int(180 - (start_angle - limits[0]) * 
                        (180.0/(limits[1]-limits[0])))
            end = int(180 - (end_angle - limits[0]) * 
                      (180.0/(limits[1]-limits[0])))
        else:
            start = int( (end_angle - limits[0]) * 
                         (180.0/(limits[1]-limits[0])))
            end = int( (start_angle - limits[0]) * 
                       (180.0/(limits[1]-limits[0])))
    
        #Draw the arch
        theta = start
        radius = 0.85
        while (theta >= end):
            x_vect.append( radius * math.cos(theta * (pi/180)) )
            y_vect.append( radius * math.sin(theta * (pi/180)) )
            theta -= 1
            
        theta = end
        radius = 1.0
        while (theta <= start):
            x_vect.append( radius * math.cos(theta * (pi/180)) )
            y_vect.append( radius * math.sin(theta * (pi/180)) )
            theta += 1
        
        if( border ):
            #Close the loop
            x_vect.append(-0.85)
            y_vect.append(0.0)
            
            p = self.plot(x_vect, y_vect, 'b-', color='black', linewidth=1.5)
        else:
            p = self.fill(x_vect, y_vect, colour, linewidth=0.0, alpha=0.4)


    def draw_needle(self, current_value, limits, graph_positive):
        x_vect = []
        y_vect = []
        
        if current_value == None:
            self.text(0.0, 0.4, "N/A", size=10, va='bottom', ha='center')
        else:
            self.text(0.0, 0.4, "%.2f" % current_value, size=10, va='bottom', 
                      ha='center')
        
            #Clamp the value to the limits
            if( current_value < limits[0] ):
                current_value = limits[0]
            if( current_value > limits[1] ):
                current_value = limits[1]
            
            theta = 0
            length = 0.95
            if( graph_positive ):
                angle = 180.0 - (current_value - limits[0])*(180.0/abs(limits[1]-limits[0]))
            else:
                angle = (current_value - limits[0])*(180.0/abs(limits[1]-limits[0]))
            
            while (theta <= 270):
                x_vect.append( length * math.cos((theta + angle) * (pi/180)) )
                y_vect.append( length * math.sin((theta + angle) * (pi/180)) )
            length = 0.05
            theta += 90
            
            p = self.fill(x_vect, y_vect, 'b', alpha=0.4)



    def draw_ticks(self, limits, graph_positive):
        if( graph_positive ):
            angle = 180.0
        else:
            angle = 0.0
        i = 0
        j = limits[0]
        
        while( i*limits[3] + limits[0] <= limits[1] ):
            x_vect = []
            y_vect = []
            if( i % (limits[2]/limits[3]) == 0 ):
                x_pos = 1.1 * math.cos( angle * (pi/180.0))
                y_pos = 1.1 * math.sin( angle * (pi/180.0))
                if( type(limits[2]) is types.FloatType ):
                    self.text( x_pos, y_pos, "%.2f" % j, size=10, va='center', 
                               ha='center', rotation=(angle - 90)) 
                else:
                    self.text( x_pos, y_pos, "%d" % int(j), size=10, 
                               va='center', ha='center', rotation=(angle - 90)) 
                tick_length = 0.15
                j += limits[2]
            else:
                tick_length = 0.05
            i += 1
            x_vect.append( 1.0 * math.cos( angle * (pi/180.0)))
            x_vect.append( (1.0 - tick_length) * math.cos( angle * (pi/180.0)))
            y_vect.append( 1.0 * math.sin( angle * (pi/180.0)))
            y_vect.append( (1.0 - tick_length) * math.sin( angle * (pi/180.0)))
            p = self.plot(x_vect, y_vect, 'b-', linewidth=1, alpha=0.4, 
                          color="black")
            if( graph_positive ):
                angle -= limits[3] * (180.0/abs(limits[1]-limits[0]))
            else:
                angle += limits[3] * (180.0/abs(limits[1]-limits[0]))
        if( i % (limits[2]/limits[3]) == 0 ):
            x_pos = 1.1 * math.cos( angle * (pi/180.0))
            y_pos = 1.1 * math.sin( angle * (pi/180.0))
            if( type(limits[2]) is types.FloatType ):
                self.text( x_pos, y_pos, "%.2f" % j, size=10, va='center', 
                           ha='center', rotation=(angle - 90)) 
            else:
                self.text( x_pos, y_pos, "%d" % int(j), size=10, va='center', 
                           ha='center', rotation=(angle - 90)) 


    def draw_bounding_box(self, graph_width, graph_height):
        x_vect = [
        graph_width/2,
        graph_width/2,
        -graph_width/2,
        -graph_width/2,
        graph_width/2,
        ]
        
        y_vect = [
        -0.1,
        graph_height,
        graph_height,
        -0.1,
        -0.1,
        ]
    
        p = self.plot(x_vect, y_vect, 'r-', linewidth=0)


if __name__=='__main__':
    from pylab import figure, show
    
    current_value = -4.0
    limits = [-1.0,1.0,1,0.1]
    zone_colour = [[-1.0,0.0,'r'],[0.0,0.5,'y'],[0.5,1.0,'g']]
    attribute_name = "Rx MOS (24h)"
    
    
    graph_height = 1.6
    graph_width = 2.4
    fig_height = graph_height
    fig_width = graph_width
    
    fig = figure(figsize=(fig_width, fig_height ))
    
    rect = [(0.0/fig_width), (0.2/fig_height),
            (graph_width/fig_width), (graph_height/fig_height)]
    
    gauge = Gauge(fig, rect,
                  xlim=( -0.1, graph_width+0.1 ),
                  ylim=( -0.4, graph_height+0.1 ),
                  xticks=[],
                  yticks=[],
                  )
    gauge.set_axis_off()
    fig.add_axes(gauge)
    
    show()
