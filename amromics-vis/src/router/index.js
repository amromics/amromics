/* eslint-disable */
import Vue from 'vue'
import Router from 'vue-router'
import SingleSample from '@/views/SingleSample'
import Collection from '@/views/Collection'
Vue.use(Router)

export default new Router({
   mode: 'history',
  routes: [{
      path: '/sample/:id',

      component: SingleSample

    },

    {
      path: '/',

      component: Collection

    }
  ],
})
