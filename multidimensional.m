clear all
metodo = 0;

while(metodo < 1 || metodo > 4)
    clear
    clc
    fprintf('Escolha o metodo de resolucao: \n');
    fprintf('1 - Gradiente \n');
    fprintf('2 - Gradiente ótimo \n');
    fprintf('3 - Newton Multidimensional \n');
    fprintf('4 - Newton Modificado \n');
    metodo = input(' ');
end

clc

syms x1 x2 alpha;
f = (x1-2)^4+(x1-2*x2)^2;

x = [x1,x2];

clc

if(metodo == 1) %E = 0.0001, alfa = 0.01, x0 = (1,1) e N = 200
   fprintf('---Método do Gradiente--- \n');
   a = input('Digite o ponto inicial x1: ');
   b = input('Digite o ponto inicial x2: ');
   eps = input('Digite o valor de tolerância: ');
   alfa = input('Digite o valor de alfa: ');
   N = input('Digite o número de iterações: ');
   k = 1;
    
   v(k,:) = [a,b];
   for i=1:2 %define a estrutura do gradiente inicial
      grad(i)=(diff(f,x(i)));
      clc;
   end
    
    gradf(k,:) = subs(grad,x,v(k,:)); %define o gradiente de f, considerando v(k,:)
    normGradf(k) = norm(gradf(k,:));  %define a norma do gradiente de f

    while(normGradf(k) >= eps && k <= N) %critério de parada do algoritmo
       d(k,:) = gradf(k,:)/normGradf(k); %define a direção de descida d(k,:)
       v(k+1,:) = v(k,:) - (alfa*d(k,:)); %atualiza o valor das variaveis a partir de alfa
       gradf(k+1,:) = subs(grad,x,v(k+1,:)); %atualiza o valor do gradiente de f considerando v(k+1,:)           
       normGradf(k+1) = norm(gradf(k+1,:)); %define a norma do gradiente atualizado
       k = k+1;
    end
    ffinal = subs(f,x,v(k,:));
    
    %Resultados
    fprintf('Pontos(x1,x2): \n');
    xf = [v(k,1) v(k,2)]
    fprintf(1, 'Valor da função = %6.4f \n', ffinal); 
    fprintf(1, 'Valor de k = %6.4f \n', k-1);
end

if(metodo == 2) %E = 0.0001, x0 = (1,1) e N = 200  
    fprintf('--Método do Gradiente Ótimo--- \n ');
    fprintf('Escolha a sub-rotina desejada: \n');
    fprintf('1 - Secao Aurea\n');
    fprintf('2 - Bissecao\n');
    fprintf('3 - Fibonacci\n');
    fprintf('4 - Newton\n');
    sub = input(' ');
    clc
    a = input('Digite o ponto inicial x1: ');
    b = input('Digite o ponto inicial x2: ');
    eps = input('Digite o valor de tolerância: ');
    N = input('Digite o número máximo de iterações: ');
    k = 1;
    
    v(k,:) = [a,b];
    for i=1:2
        grad(i)=(diff(f,x(i))); %define a estrutura de gradiente de f
        clc;
    end
    
    gradf(k,:) = subs(grad,x,v(k,:)); %atualiza o gradiente de f a partir de v(k,:)
    normGradf(k) = norm(gradf(k,:));  %define a norma do gradiente atualizado
    
    while(normGradf(k) >= eps && k <= N) %critério de parada do algoritmo
        d(k,:) = gradf(k,:)/normGradf(k); %define a direção d(k,:)
        vb(k,:) = v(k,:)-(alpha*d(k,:)); %guarda as variáveis atualizadas em relação à v(k,:)-alpha*d(k,:)
        fb(k) = subs(f,x,vb(k,:)); %guarda o valor de f em função de vb(k,:) - com sym alpha
        kn = 1;
        
        if(sub == 1) 
            %Secao aurea - Inicio
            a = 0;
            b = 1;
            lambda = a + (1-0.618) * (b-a); %atribuindo valor de lambda
            mi = a + 0.618 * (b-a); %atribuindo valor de mi
    
            %Passo principal
            while(b-a >= eps) %critério de parada do método
               if(subs(fb,alpha,lambda) > subs(fb,alpha,mi)) 
                  a = lambda;
                  lambda = mi;
                  mi = a + 0.618*(b-a);
               else 
                  b = mi;
                  mi = lambda;
                  lambda = a + (1-0.618) * (b-a);
               end
            end
            alfa = (a+b)/2;
        end
        %Secao Aurea - Fim
        
        %Bissecao - Inicio
        if(sub == 2)    
           %Dados Iniciais
           fb_d1 = diff(fb(k)); %derivada primeira de fb(k)
           n = log(eps/(b(1)-a(1)))/log(1/2); %cálculo de n
           a = 0;
           b = 1;
           %Passo Principal
           while 1
                lambda = 1/2*(b+a); %calculando valor de lambda
                teta_lambda_d1 = subs(fb_d1,alpha,lambda); %calculo da derivada primeira de teta lambda
                fd1 = eval(teta_lambda_d1);
                if(fd1 == 0 || kn > n) %condição de parada do método
                    break;
                elseif(fd1 > 0)               
                    b = lambda;
                elseif(fd1 < 0)
                    a = lambda; 
                end
                kn = kn+1;
           end
           alfa = (a+b)/2; %valor apropriado para alfa
        end

        %bissecao - fim
        

        %fibonacci - inicio
        if(sub == 3)

            %Passo inicial

            n = 0;
            D = 0.001;
            af(kn) = 0;
            bf(kn) = 1;
            
            while(fibonacci(n) < (bf(kn)-af(kn))/eps) %define o tamanho de n
                n = n+1;
            end
            
            lambda(kn) = af(kn) + fibonacci(n-2)*(bf(kn)-af(kn))/fibonacci(n); %define o valor inicial de lambda
            mi(kn) = af(kn) + fibonacci(n-1)*(bf(kn)-af(kn))/fibonacci(n); %define o valor inicial de mi
 
            %Passo Principal
            while 1
                if(subs(fb,alpha,lambda(kn)) > subs(fb,alpha,mi(kn))) 
                    af(kn+1) = lambda(kn);
                    bf(kn+1) = bf(kn);
                    lambda(kn+1) = mi(kn);
                    mi(kn+1) = af(kn+1) + fibonacci(n-kn-1)*(bf(kn+1)-af(kn+1))/fibonacci(n-kn); %atualiza valor de mi utilizando a sequência de fibonacci
                else
                    af(kn+1) = af(kn);
                    bf(kn+1) = mi(kn);
                    mi(kn+1) = lambda(kn);
                    lambda(kn+1) = af(kn+1) + fibonacci(n-kn-2)*(bf(kn+1)-af(kn+1))/fibonacci(n-kn); %atualiza valor de lambda utilizando a sequêndia de fibonacci
                end
                if(kn == n-2) %critério de parada do método
                    lambda(n) = lambda(n-1);
                    mi(n) = lambda(n-1)+D;
                    if(subs(fb,alpha,lambda(n)) > subs(fb,alpha,mi(n)))
                        af(n) = lambda(n);
                        bf(n) = bf(n-1);
                    else
                        af(n) = af(n-1);
                        bf(n) = lambda(n);
                    end
                    break;
                end
                kn = kn+1;
            end       
            alfa = (af(n)+bf(n))/2; %valor apropriado para alfa
        end
        %fibonacci - fim
        
        %newton - inicio
        if(sub == 4)
           fb_d1 = diff(fb(k)); %derivada primeira de fb(k)
           fb_d2 = diff(fb(k),2); %derivada segunda de fb(k) 
           lambda(1) = 6;
            while 1
                teta_lambda_d1 = subs(fb_d1,alpha,lambda(kn));
                teta_lambda_d2 = subs(fb_d2,alpha,lambda(kn));
                fd1 = eval(teta_lambda_d1);
                fd2 = eval(teta_lambda_d2);
                lambda(kn+1) = lambda(kn) - (fd1/fd2);
            
                if(abs(lambda(kn+1) - lambda(kn)) < eps)
                  break;
                else
                  kn = kn+1;
                end
            end
            alfa = lambda(kn);
        end
        %newton - fim
        
        v(k+1,:) = v(k,:) - alfa*d(k,:); %atualiza as variáveis em função de lambda(kn), encontrada pelo método de newton
        ffinal = subs(f,x,v(k+1,:)); %atualiza o valor de saída da função objetivo
        gradf(k+1,:) = subs(grad,x,v(k+1,:)); %atualiza o gradiente de f
        normGradf(k+1,:) = norm(gradf(k+1,:)); %guarda a norma do gradiente atualizado de f
        k = k+1;
    end
    
    %Resultados
    fprintf('Pontos (x1,x2): \n');
    xf = [v(k,1) v(k,2)]
    fprintf(1, 'Valor da Função = %6.4f \n', ffinal);
    fprintf(1, 'Valor de k = %6.4f \n', k-1);
end

if(metodo == 3) %E = 0.0001, x0 = (1,1) e N = 200
    fprintf('Método de Newton Multidimensional \n');
    a = input('Digite o ponto inicial x1: ');
    b = input('Digite o ponto inicial x2: ');
    eps = input('Digite um valor de tolerância positivo: ');
    N = input('Digite um número máximo de iterações N: ');
    k = 1;
    
    %vetor de variáveis xk
    v(k,:) = [a,b];
    
    %define a estrutura do gradiente de f
    for m = 1:2
        gradf(m) = diff(f,x(m));
    end
    
    %define a estrutura da matriz hessiana de f
    for i = 1:2
        for j = 1:2
            H(i,j) = diff(gradf(i),x(j));
        end
    end
   
    %Passo princial
    while(k <= N)
        ffinal(k) = subs(f,x,v(k,:)); %atualiza o resultado da função objetivo
        
        %vetor gradiente
        for p = 1:2
            vgradf(k,p) = subs(gradf(p),x,v(k,:));
        end
        
        %calcula o módulo do vetor gradiente
        normGradf(k) = norm(vgradf(k,:));
        
        %atualiza a matriz hessiana
        h(:,:,k) = subs(H,x,v(k,:));
        
        %calcula a inversa da matriz hessiana atualizada
        invH(:,:,k) = inv(h(:,:,k));
        
        %xk+1 = xk - H(xk)^-1 * gradf(xk)
        v(k+1,:) = v(k,:) - (invH(:,:,k)*vgradf(k,:)')'; 
        
        if(normGradf(k) < eps)
            break;
        else
             k = k+1;
        end
    end
    
    %Resultados
    fprintf('Pontos [x1,x2]: \n');
    v(k,:)
    fprintf(1, 'Valor da Função = %6.4f \n', ffinal(k));
    fprintf(1, 'Valor de k = %6.4f \n', k-1);
end

if(metodo == 4) %x0 = (1,1), N = 200 e E pequeno = 0.0001
    fprintf('---Método Hibrido--- \n');
    a = input('Digite o ponto inicial x1: ');
    b = input('Digite o ponto inicial x2: ');
    eps = input('Digite o valor de tolerância: ');
    N = input('Digite o número máximo de iterações: ');
    k = 1;
    
    %vetor de variáveis xk
    v(k,:) = [a,b];
    
    %define a estrutura do gradiente de f
    for m = 1:2
        gradf(m) = diff(f,x(m));
    end
    
    %define a estrutura da matriz hessiana de f
    for i = 1:2
        for j = 1:2
            H(i,j) = diff(gradf(i),x(j));
        end
    end
    
    %Passo princial
    while(k <= N)
        if(k >= 4)
            eps2 = eps;
        else
            eps2 = 100;
        end
        
        %vetor gradiente
         for p = 1:2
           vgradf(k,p) = subs(gradf(p),x,v(k,:));
         end
        
         %calcula o módulo do vetor gradiente
         normGradf(k) = norm(vgradf(k,:));
        
         %atualiza a matriz hessiana
         h(:,:,k) = subs(H,x,v(k,:));
        
         %calcula a inversa da matriz hessiana atualizada
         invH(:,:,k) = inv(h(:,:,k));
        
         %M(k) = [E(k)*I+H(xk)]^-1
         mm(:,:,k) = eps2*invH(:,:,k) + h(:,:,k);
         M(:,:,k) = inv(mm(:,:,k));

         %xk+1 = xk - M(k) * gradf(xk)
         v(k+1,:) = v(k,:) - (M(:,:,k)*vgradf(k,:)')'; 
        
         if(normGradf(k) < eps)
            break;
         else
            k = k+1;
         end
    end
    
    %Resultados
    fprintf('Pontos [x1,x2]: \n');
    v(k,:)
    fprintf(1, 'Valor da Função = %6.4f \n', subs(f,x,v(k,:))); %resultado da função objetivo);
    fprintf(1, 'Valor de k = %6.4f \n', k-1);
end
    


