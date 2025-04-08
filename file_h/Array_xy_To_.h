#pragma once

//template<typename T>sf::VertexArray Array_xy_To_(std::vector<std::pair<T, T>>& Array_xy) {
//	sf::VertexArray 
//
//	for (auto const& : Array_xy) {
//        for (float x = 0; x < 10.0f; x += 0.1f) {
//            float y = F((float)x);
//            float screenX = x * xScale;
//            float screenY = yOffset - y * amplitude;
//            graph.append(sf::Vertex{ sf::Vector2f(screenX, screenY),sf::Color::Red });
//        }
//	}
//
//}

#if __has_include(<SFML/Graphics.hpp>)
#include <SFML/Graphics.hpp>
#include <cmath>
#include <functional>

const sf::Vector2u
Drawing_Window_Size = { 800,600 },
center_Window = { Drawing_Window_Size.x / 2 -200 ,Drawing_Window_Size.y / 2 + 200 };
const sf::Vector2f ratioX_to_Y = { 20.0f,12.5f };
const double modifire = 1;


//
//funct(F&& f,
// const sf::Vector2u Drawing_Window_Size = { 800,600 },
// const sf::Vector2u center_Window = { Drawing_Window_Size.x / 2 - 200 ,Drawing_Window_Size.y / 2 + 200 },
// const sf::Vector2f ratioX_to_Y = { 20.0f,12.5f },
// const double modifire = 1)
template<typename F>sf::VertexArray funct(F&& f,
const sf::Vector2u Drawing_Window_Size = { 800,600 },
const sf::Vector2u center_Window = { Drawing_Window_Size.x / 2 - 200 ,Drawing_Window_Size.y / 2 + 200 },
const sf::Vector2f ratioX_to_Y = { 20.0f,12.5f },
const double modifire = 1) {
    sf::VertexArray graph(sf::PrimitiveType::LineStrip, 0);

    const float xScale = ratioX_to_Y.x * modifire;
    const float yScale = ratioX_to_Y.y * modifire;
    const float yOffset = center_Window.y;

    for (float x = -10.f; x < 10.0f; x += 0.1f) {
        auto y = f(static_cast<double>(x));  
        const float screenX = x * xScale + center_Window.x;
        const float screenY = yOffset - static_cast<float>(y) * yScale;
        graph.append(sf::Vertex{ sf::Vector2f(screenX, screenY), sf::Color::Red });
    }
    return graph;
}

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
        window.draw(funct((f)));
        window.display();
    }

}
#else
    #pragma message("SFML not found: Skipping related code")
#endif