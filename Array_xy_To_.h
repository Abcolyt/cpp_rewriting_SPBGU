#pragma once
//The file of functions for rendering the function y=f(x)
//#define NOMINMAX
//#define WIN32_LEAN_AND_MEAN
//#include <windows.h>
//#undef ERROR
//#undef IGNORE
//#undef IN
//#undef OUT

#if __has_include(<SFML/Graphics.hpp>)
#include <SFML/Graphics.hpp>  
#include <SFML/Window/Keyboard.hpp>
#include <SFML/Graphics.hpp>

#include <vector>
#include <functional>
#include <cmath>
#include <string>

namespace Drawing_const {};
using namespace Drawing_const;
namespace Drawing_const {
    
    //window size
    sf::Vector2u Drawing_Window_Size = { 1600,900 };
    //the center of the screen (the center of the graph, that is, the default point (0,0))
    sf::Vector2u center_Window = { Drawing_Window_Size.x / 2,Drawing_Window_Size.y / 2 };

    //the ratio of relative x and y coordinates
    sf::Vector2f ratio_of_modifiers_by_xy = { 16.0f,9.0f };

    //the vector of the graph offset from the center of the screen
    sf::Vector2f  offset_from_Zero_vector = { 0,0 };
    //scale
    float modifire = 1;

    float pointRadius = 3.0f; // R of point
    float risk_lenght = 20.0;

    double left_border  = (-1) * ((center_Window.x) / (ratio_of_modifiers_by_xy.x * modifire)),
           right_border = (Drawing_Window_Size.x - center_Window.x) / (ratio_of_modifiers_by_xy.x * modifire),
           step_by_x = 0.1 / modifire;
    void update_borders() {
        left_border = (-1) * ((center_Window.x) / (ratio_of_modifiers_by_xy.x * modifire));
        right_border = (Drawing_Window_Size.x - center_Window.x) / (ratio_of_modifiers_by_xy.x * modifire);
        step_by_x = 0.1 /( modifire * ratio_of_modifiers_by_xy.x);
    }

}


//
//funct(F&& f,,double a=0, double b=10,
// const sf::Vector2u Drawing_Window_Size = { 1600,900 },
// const sf::Vector2u center_Window = { Drawing_Window_Size.x / 2 - 200 ,Drawing_Window_Size.y / 2 + 200 },
// const sf::Vector2f ratioX_to_Y =  { 16.0f,9.0f },
// const double modifire = 1)
template<typename F>sf::VertexArray funct(F&& f,double a= left_border, double b= right_border
/*,const sf::Vector2u Drawing_Window_Size = { 800,600 },
const sf::Vector2u center_Window = { Drawing_Window_Size.x / 2  ,Drawing_Window_Size.y / 2 },
const sf::Vector2f ratioX_to_Y = { 16.0f,9f },
const double modifire = 1*/) {
    //std::cout << "b="<<b;
    sf::VertexArray graph(sf::PrimitiveType::LineStrip, 0);

    const float xScale = ratio_of_modifiers_by_xy.x * modifire;
    const float yScale = ratio_of_modifiers_by_xy.y * modifire;

    
    for (float x = a; x <= b; x += step_by_x) {
        auto y = f(static_cast<double>(x- offset_from_Zero_vector.x))- offset_from_Zero_vector.y;
        const float screenX = x * xScale + center_Window.x;
        const float screenY = center_Window.y - static_cast<float>(y) * yScale;
        graph.append(sf::Vertex{ sf::Vector2f(screenX, screenY), sf::Color::Red });
    }
    return graph;
}



template<class T>void DrawFunctions(const std::vector<std::function<T(T)>>& functions, T a = (T)0, T b = (T)0, const std::vector<std::pair<T, T>>& points = std::vector<std::pair<T, T>>{}, const std::string& window_name = "y(x)=", const std::vector<std::function<T(T)>>& cyan_functions = std::vector<std::function<T(T)>>{}) {
    sf::RenderWindow window(sf::VideoMode(Drawing_Window_Size), window_name);

    window.setVerticalSyncEnabled(true);
    //sf::ContextSettings settings = window.getSettings();
    //glPointSize(5);



    while (window.isOpen()) {

        // // // // make coordinate axes
        sf::VertexArray axes(sf::PrimitiveType::Lines, 4);

        //Calculating the scales along the X and Y axes (as in the function)
        const float xScale = ratio_of_modifiers_by_xy.x * modifire;
        const float yScale = ratio_of_modifiers_by_xy.y * modifire;

        //  OX (y = 0 in relative coordinates)
        float screenY_OX = center_Window.y - (0 - offset_from_Zero_vector.y) * yScale;
        axes.append(sf::Vertex(sf::Vector2f(0, screenY_OX), sf::Color::Black));
        axes.append(sf::Vertex(sf::Vector2f(Drawing_Window_Size.x, screenY_OX), sf::Color::Black));

        //  OY (x = 0 in relative coordinates)
        float screenX_OY = center_Window.x - (0 - offset_from_Zero_vector.x) * xScale;
        axes.append(sf::Vertex(sf::Vector2f(screenX_OY, 0), sf::Color::Black));
        axes.append(sf::Vertex(sf::Vector2f(screenX_OY, Drawing_Window_Size.y), sf::Color::Black));
        //
        // // // //
        
        // // // //createborder markings
       
        sf::VertexArray markers(sf::PrimitiveType::Lines, 4); 

        float screenX_a = center_Window.x + (static_cast<float>(a) + offset_from_Zero_vector.x) * xScale;
        const float risk_size = risk_lenght;
        markers.append(sf::Vertex(sf::Vector2f(screenX_a, screenY_OX - risk_size/2), sf::Color::Black));
        markers.append(sf::Vertex(sf::Vector2f(screenX_a, screenY_OX + risk_size/2), sf::Color::Black));

        float screenX_b = center_Window.x + (static_cast<float>(b) + offset_from_Zero_vector.x) * xScale;
        markers.append(sf::Vertex(sf::Vector2f(screenX_b, screenY_OX - risk_size/2), sf::Color::Black));
        markers.append(sf::Vertex(sf::Vector2f(screenX_b, screenY_OX + risk_size/2), sf::Color::Black));
        // // // //


        // // // //make function lines
        std::vector<sf::VertexArray> graphs;
        for (size_t i = 0; i < functions.size(); ++i) {
            graphs.push_back(funct(functions[i])); 

            if (i == 0) {
                for (size_t j = 0; j < graphs.back().getVertexCount(); ++j) {
                    graphs.back()[j].color = sf::Color::Black;
                }
            }

            if ((i == functions.size() - 1)&&(functions.size()!=1)) {
                for (size_t j = 0; j < graphs.back().getVertexCount(); ++j) {
                    graphs.back()[j].color = sf::Color::Blue;
                }
            }
        }
        if (cyan_functions.size() != 0) {
            for (size_t i = 0; i < cyan_functions.size(); ++i) {
                graphs.push_back(funct(cyan_functions[i]));
                for (size_t j = 0; j < graphs.back().getVertexCount(); ++j) {
                    graphs.back()[j].color = sf::Color::Cyan;
                }
            }
        }
        // // // //

        // // // // point array
        std::vector<sf::CircleShape> pointShapes;
        
        for (const auto& pt : points) {
            // mathematic -> screen coordinate
            float screenX = center_Window.x + (static_cast<float>(pt.first) + offset_from_Zero_vector.x) * xScale;
            float screenY = center_Window.y - (static_cast<float>(pt.second) - offset_from_Zero_vector.y) * yScale;

            sf::CircleShape point(pointRadius);
            point.setFillColor(sf::Color::Green);

            point.setOrigin(sf::Vector2f(pointRadius, pointRadius));
            point.setPosition(sf::Vector2f(screenX, screenY));

            pointShapes.push_back(point);
        }
        // // // //

        // // // // Event handling
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>() || (event->is<sf::Event::KeyPressed>() && event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::Escape))
            {
                window.close();
            }

            if (event->is<sf::Event::KeyPressed>()) {
#define MAKRO_SHIFT(a,b) Drawing_const::offset_from_Zero_vector = Drawing_const::offset_from_Zero_vector +  sf::Vector2f{ a/modifire,b/modifire };Drawing_const::update_borders();
#define MAKRO_ZOOM(k) modifire = modifire * k;Drawing_const::step_by_x*=k;Drawing_const::update_borders(); // pointRadius*=k;
#define RATIO_MODIFIRE(mx,my)   ratio_of_modifiers_by_xy = sf::Vector2f{mx* ratio_of_modifiers_by_xy.x ,ratio_of_modifiers_by_xy.y*my};Drawing_const::update_borders();
                
                using namespace sf::Keyboard;
                switch (event->getIf<sf::Event::KeyPressed>()->code)
                {
                case Key::Escape:window.close(); break;
                case Key::A:MAKRO_SHIFT(16 / ratio_of_modifiers_by_xy.x, 0); break;
                case Key::D:MAKRO_SHIFT(-16 / ratio_of_modifiers_by_xy.x, 0); break;
                case Key::W:MAKRO_SHIFT(0, 9 / ratio_of_modifiers_by_xy.y); break;
                case Key::S:MAKRO_SHIFT(0, -9 / ratio_of_modifiers_by_xy.y); break;
                case Key::Space: MAKRO_ZOOM(1.5); break;
                case Key::LAlt: MAKRO_ZOOM(0.333333); break;
                case Key::Left:   RATIO_MODIFIRE(2.f / 3.f, 1); break;
                case Key::Right:  RATIO_MODIFIRE(3.f / 2.f, 1); break;
                case Key::Down:   RATIO_MODIFIRE(1, 2.f / 3.f); break;
                case Key::Up:     RATIO_MODIFIRE(1, 3.f / 2.f); break;

                case Key::Tab:
                    std::cout << "current rendering parameters:\n"
                        << "l:" << left_border << "\n"
                        << "r:" << right_border << "\n"
                        << "m:" << modifire << "\n"
                        << "s:" << step_by_x << ";\n"
                        << "camera:(" << -offset_from_Zero_vector.x << ";" << offset_from_Zero_vector.y << ")\n"
                        << "ratio_:{" << ratio_of_modifiers_by_xy.x << ":" << ratio_of_modifiers_by_xy.y << "}\n"
                        << "pointRadius:{" << pointRadius << "}\n"; break;

                case Key::F1:
                    std::cout << "rendering keys\n"
                        << "A:MAKRO_SHIFT(1,0)\n"
                        << "D:MAKRO_SHIFT(-1, 0)\n"
                        << "W:MAKRO_SHIFT(0,1)\n"
                        << "S:MAKRO_SHIFT(0,-1)\n"
                        << "Space: MAKRO_ZOOM(3/2) \n"
                        << "LAlt: MAKRO_ZOOM(2 / 3) \n"
                        << "Left:graph compression by x\n"
                        << "Right:graph decompression by x\n"
                        << "Down:graph compression by y\n"
                        << "Up:graph decompression by y\n"
                        << "Tab:current rendering parameters (l:left_border ,r:right_border,m:modifire)\n"; break;
                default:std::cout << "error: there is no fixed behavior use f1 for help\n";
                }
            }
        }
        // // // // 

        // // // //  drawing
        // // // clear
        window.clear(sf::Color::White);
        // // //
 

        // // // functions
        for (const auto& graph : graphs) {
            window.draw(graph);
        }
        // // //

        // // // dots
        for (const auto& point : pointShapes) {
            window.draw(point);
        }
        // // //

        // // // lines
        window.draw(axes);
        window.draw(markers); 
        // // //

        // // // displaing
        window.display();
        // // //
        
        // // // // 

    }
}


////vector fun object
//template<class T>void DrawFunctions(const std::vector<std::function<T(T)>>& functions, const std::string& window_name = "y(x)=") {
//    DrawFunctions(functions,  std::vector<std::pair<T, T>>{}, (T)0, (T)0, window_name);
//}
#if 0
//ones fun object
template<class T>void DrawFunction(T& f) {
    sf::RenderWindow window(sf::VideoMode(Drawing_Window_Size), "y(x)=");

    // Создаём оси координат один раз перед циклом
    sf::VertexArray axes(sf::PrimitiveType::Lines, 4);

    // Ось X (горизонтальная линия)
    axes.append({ sf::Vertex{ sf::Vector2f(0         , center_Window.y), sf::Color::Black } });
    axes.append({ sf::Vertex{ sf::Vector2f(Drawing_Window_Size.x, center_Window.y), sf::Color::Black } });

    // Ось Y (вертикальная линия)
    axes.append({ sf::Vertex{ sf::Vector2f(center_Window.x ,          0), sf::Color::Black } });
    axes.append({ sf::Vertex{ sf::Vector2f(center_Window.x , Drawing_Window_Size.y), sf::Color::Black } });

    auto drawning_function = funct(f);
    while (window.isOpen()) {

        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>() ||
                (event->is<sf::Event::KeyPressed>() &&
                    event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::Escape))
            {
                window.close();
            }
        }

        window.clear(sf::Color::White);

        window.draw(axes);
        window.draw(drawning_function);
        window.display();
    }

}
#endif

template<class T, typename F>
void DrawFunction(F&& func) {
    DrawFunctions<T>({ std::function<T(T)>(std::forward<F>(func)) });
}

#else
    #pragma message("SFML not found: Skipping related code")
#endif