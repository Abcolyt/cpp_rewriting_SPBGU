#pragma once



#if __has_include(<SFML/Graphics.hpp>)
#include <SFML/Graphics.hpp>
#include <vector>
#include <functional>
#include <cmath>

namespace Drawing_const {};
using namespace Drawing_const;
namespace Drawing_const {
    sf::Vector2u
        Drawing_Window_Size = { 1600,900 },
        center_Window = { Drawing_Window_Size.x / 2,Drawing_Window_Size.y / 2 };
    
    sf::Vector2f ratio_of_modifiers_by_xy = { 16.0f,9.0f },
        offset_from_Zero_vector = { 0,0 };
    const double modifire = 1;

    const double left_border = (-1) * ((center_Window.x) / (ratio_of_modifiers_by_xy.x * modifire)),
        right_border = (Drawing_Window_Size.x - center_Window.x) / (ratio_of_modifiers_by_xy.x * modifire),
        step_by_x = 0.1;
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

//ones fun object
template<class T>void draw_function(T f) {
    sf::RenderWindow window(sf::VideoMode(Drawing_Window_Size), "y(x)=");
    
    // Создаём оси координат один раз перед циклом
    sf::VertexArray axes(sf::PrimitiveType::Lines, 4);

    // Ось X (горизонтальная линия)
    axes.append({ sf::Vertex{ sf::Vector2f(         0         , center_Window.y), sf::Color::Black } });
    axes.append({ sf::Vertex{ sf::Vector2f(Drawing_Window_Size.x, center_Window.y), sf::Color::Black } });

    // Ось Y (вертикальная линия)
    axes.append({ sf::Vertex{ sf::Vector2f(center_Window.x ,          0         ), sf::Color::Black } });
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

//vector fun object
template<class T>void draw_functions(const const std::vector<std::function<T(T)>>& functions) {
    sf::RenderWindow window(sf::VideoMode(Drawing_Window_Size), "y(x)=");
    sf::VertexArray axes(sf::PrimitiveType::Lines, 4);

    // OX 
    axes.append({ sf::Vertex{ sf::Vector2f(0         , center_Window.y), sf::Color::Black } });
    axes.append({ sf::Vertex{ sf::Vector2f(Drawing_Window_Size.x, center_Window.y), sf::Color::Black } });
    // OY 
    axes.append({ sf::Vertex{ sf::Vector2f(center_Window.x ,          0), sf::Color::Black } });
    axes.append({ sf::Vertex{ sf::Vector2f(center_Window.x , Drawing_Window_Size.y), sf::Color::Black } });

    
    std::vector<sf::VertexArray> graphs;
    for (const auto& func : functions) {
        graphs.push_back(funct(func));
    }
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
        for (const auto& graph : graphs) {
            window.draw(graph);
        }
        window.display();
    }
}


#else
    #pragma message("SFML not found: Skipping related code")
#endif